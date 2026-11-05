# =============================================================================
# 02b Genotypes of one chromosome for the tissue's donors (from the plink file of 02a)
#
# Variants of the chromosome with MAF >= MAF_THRESHOLD among the tissue's
# donors, read from the plink .bed and stored as raw codes (0, 1, 2 copies of
# the second .bim allele, the VCF REF allele for a plink2 conversion; 3 =
# missing). Independent of 01 (the donor list is recomputed from the same
# inputs), so 01 and 02a/02b can run at the same time.
#
# Run:     sbatch 02b_prepareGenotypes.sh        (array over chromosomes 1-22)
# Usage:   Rscript 02b_prepareGenotypes.R <chromosome>
# Outputs: GENO_DIR/<tissue>/genotypes_chr<c>.rds: G (donors x variants, raw),
#          variants (id, pos, a1, a2, maf)
# =============================================================================

.a <- commandArgs(FALSE); .f <- sub('^--file=', '', .a[grep('^--file=', .a)])
source(file.path(if (length(.f)) dirname(normalizePath(.f)) else '.', '00_config.R'))
chrom <- as.integer(commandArgs(trailingOnly = TRUE)[1])

donors <- tissueDonors(TISSUE)
fam <- fread(paste0(GENOTYPE_PLINK, '.fam'), header = FALSE)
keep <- match(donors, gtexDonorId(fam$V2))
bim <- readBim(GENOTYPE_PLINK)[, row := .I][sub('^chr', '', chr) == as.character(chrom)]
message('chr', chrom, ': ', nrow(bim), ' variants, ', length(donors), ' donors')
if (!nrow(bim)) stop('No variants on chromosome ', chrom)

G <- readPlinkChromosome(GENOTYPE_PLINK, bim$row, nrow(fam), keep)
codes <- matrix(as.integer(G), nrow(G))
codes[codes == 3] <- NA
af <- colMeans(codes, na.rm = TRUE) / 2
rm(codes)
ok <- !is.na(af) & pmin(af, 1 - af) >= MAF_THRESHOLD
G <- G[, ok, drop = FALSE]
variants <- bim[ok, .(id, pos, a1, a2, maf = pmin(af, 1 - af)[ok])]
o <- order(variants$pos)
G <- G[, o, drop = FALSE]; variants <- variants[o]
dimnames(G) <- list(donors, variants$id)
message(ncol(G), ' variants with MAF >= ', MAF_THRESHOLD)
saveRDS(list(G = G, variants = variants), GENO_FILE(chrom), compress = FALSE)
