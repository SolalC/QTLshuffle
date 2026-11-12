# =============================================================================
# 02a GTEx WGS VCF to plink, once per tissue
#
#   1. One pass over GENOTYPE_VCF with plink2, keeping the tissue's donors (those
#      with expression, covariates and a DIP), FILTER = PASS, biallelic variants
#      on chr1-22, MAF >= MAF_THRESHOLD among those donors; half-calls set to
#      missing; every variant renamed chr:pos:REF:ALT (alleles truncated at
#      50 bp), since many VCF IDs are missing ('.'); variants whose new ID is
#      duplicated dropped (later steps index variants by ID). Skipped if the
#      plink file already exists.
#   2. The plink file split by chromosome (plink2 --chr), for plink2 --clump
#      in 04 (a whole-genome .bim per clumping run would be slow).
# 03-05 read the whole-genome file with BEDMatrix.
#
# Run:     sbatch 02a_vcfToPlink.sh
# Usage:   Rscript 02a_vcfToPlink.R <threads> <memory MB>
# Outputs: GENOTYPE_PLINK.{bed,bim,fam,log,keep}; PLINK_CHR(1..22).{bed,bim,fam,log}
# =============================================================================

.a <- commandArgs(FALSE); .f <- sub('^--file=', '', .a[grep('^--file=', .a)])
source(file.path(if (length(.f)) dirname(normalizePath(.f)) else '.', '00_config.R'))
args <- commandArgs(trailingOnly = TRUE)
threads <- args[1]; memory <- args[2]

plink <- function(...) {
  status <- system2(PLINK2, c(..., '--threads', threads, '--memory', memory))
  if (status != 0) stop('plink2 failed: ', paste(c(...), collapse = ' '))
}

if (file.exists(paste0(GENOTYPE_PLINK, '.fam'))) {
  message(GENOTYPE_PLINK, '.fam exists: VCF conversion skipped')
} else {
  dir.create(dirname(GENOTYPE_PLINK), recursive = TRUE, showWarnings = FALSE)
  samples <- genotypeSampleIds()             # VCF header: the .fam does not exist yet
  donors  <- tissueDonors(TISSUE)
  keep <- samples[gtexDonorId(samples) %in% donors]
  message(TISSUE, ': ', length(keep), ' of ', length(samples), ' VCF samples kept')
  keepFile <- paste0(GENOTYPE_PLINK, '.keep')
  writeLines(c('#IID', keep), keepFile)
  plink('--vcf', GENOTYPE_VCF, '--vcf-half-call', 'missing',
        '--keep', keepFile, '--var-filter', '--max-alleles', '2',
        '--chr', '1-22', '--maf', MAF_THRESHOLD,
        '--set-all-var-ids', shQuote('chr@:#:$r:$a'),  # quoted: system2 runs a shell
        '--new-id-max-allele-len', '50', 'truncate',
        '--rm-dup', 'exclude-all',
        '--make-bed', '--out', GENOTYPE_PLINK)
}

dir.create(dirname(PLINK_CHR(1)), recursive = TRUE, showWarnings = FALSE)
for (chrom in CHROMOSOMES) plink('--bfile', GENOTYPE_PLINK, '--chr', chrom, '--make-bed', '--out', PLINK_CHR(chrom))
message(nrow(fread(paste0(GENOTYPE_PLINK, '.bim'), header = FALSE)), ' variants; split into ',
        length(CHROMOSOMES), ' chromosomes')
