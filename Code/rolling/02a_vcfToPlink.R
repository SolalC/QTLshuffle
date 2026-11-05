# =============================================================================
# 02a GTEx WGS VCF to plink, once per tissue
#
# One pass over GENOTYPE_VCF with plink2, keeping the tissue's donors (those
# with expression, covariates and a DIP), FILTER = PASS, biallelic variants on
# chr1-22, MAF >= MAF_THRESHOLD among those donors; half-calls set to missing;
# every variant renamed chr:pos:REF:ALT (alleles truncated at 50 bp), since
# many VCF IDs are missing ('.'), and variants whose new ID is duplicated
# dropped (later steps index variants by ID).
# Writes a plink 1 .bed/.bim/.fam, in VCF order (sorted by chromosome and
# position), that 02b reads one chromosome at a time.
#
# Run:     sbatch 02a_vcfToPlink.sh
# Usage:   Rscript 02a_vcfToPlink.R <plink2 executable> <threads> <memory MB>
# Outputs: GENOTYPE_PLINK.{bed,bim,fam,log}; GENOTYPE_PLINK.keep (sample IDs)
# =============================================================================

.a <- commandArgs(FALSE); .f <- sub('^--file=', '', .a[grep('^--file=', .a)])
source(file.path(if (length(.f)) dirname(normalizePath(.f)) else '.', '00_config.R'))
args <- commandArgs(trailingOnly = TRUE)
plink2 <- args[1]; threads <- args[2]; memory <- args[3]

if (file.exists(paste0(GENOTYPE_PLINK, '.fam'))) stop(GENOTYPE_PLINK, '.fam exists; delete it to rebuild')
dir.create(dirname(GENOTYPE_PLINK), recursive = TRUE, showWarnings = FALSE)
samples <- genotypeSampleIds()             # VCF header: the .fam does not exist yet
donors  <- tissueDonors(TISSUE)
keep <- samples[gtexDonorId(samples) %in% donors]
message(TISSUE, ': ', length(keep), ' of ', length(samples), ' VCF samples kept')
keepFile <- paste0(GENOTYPE_PLINK, '.keep')
writeLines(c('#IID', keep), keepFile)

status <- system2(plink2, c('--vcf', GENOTYPE_VCF, '--vcf-half-call', 'missing',
                            '--keep', keepFile, '--var-filter', '--max-alleles', '2',
                            '--chr', '1-22', '--maf', MAF_THRESHOLD, 
                            '--set-all-var-ids', shQuote('chr@:#:$r:$a'),  # quoted: system2 runs a shell
                            '--new-id-max-allele-len', '50', 'truncate',
                            '--rm-dup', 'exclude-all',
                            '--make-bed', '--out', GENOTYPE_PLINK,
                            '--threads', threads, '--memory', memory))
if (status != 0) stop('plink2 failed; see ', GENOTYPE_PLINK, '.log')
message(nrow(fread(paste0(GENOTYPE_PLINK, '.bim'), header = FALSE)), ' variants written')
