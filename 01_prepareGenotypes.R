# =============================================================================
# 01  Genotype preparation (one job)
#
# Reads the GTEx WGS plink set once and writes one dosage matrix per autosome
# (donors x variants, sorted by position) for the donors with a DIP, keeping
# variants with MAF >= MAF_THRESHOLD across those donors. Tissue-level MAF and
# within-window minor allele counts are applied later (02, 03).
#
# Run:     sbatch 01_prepareGenotypes.sh
# Usage:   Rscript 01_prepareGenotypes.R
# Outputs: Results/genotypes/chr<1-22>.rds  list(G, info)
# =============================================================================

.a <- commandArgs(FALSE); .f <- sub('^--file=', '', .a[grep('^--file=', .a)])
source(file.path(if (length(.f)) dirname(normalizePath(.f)) else '.', '00_config.R'))

dip   <- loadDIP()
plink <- loadGenotypes()
X <- plink$X
colnames(X) <- gtexDonorId(colnames(X))
bim <- as.data.table(plink$bim)[, .(chr = as.character(chr), id, pos)]
rm(plink); gc()

bim[, chr := ifelse(startsWith(chr, 'chr'), chr, paste0('chr', chr))]
donors <- intersect(colnames(X), dip$SUBJID)
keep   <- !duplicated(bim$id) | bim$id == '.'
message(length(donors), ' donors with genotypes and a DIP')

for (chr in paste0('chr', 1:22)) {
  idx <- which(bim$chr == chr & keep)
  if (!length(idx)) next
  idx <- idx[order(bim$pos[idx])]
  G   <- t(X[idx, donors, drop = FALSE])
  maf <- colMeans(G, na.rm = TRUE) / 2
  ok  <- pmin(maf, 1 - maf) >= MAF_THRESHOLD
  ok[is.na(ok)] <- FALSE
  saveRDS(list(G = G[, ok, drop = FALSE], info = bim[idx][ok]),
          paste0(GENO_DIR, chr, '.rds'), compress = FALSE)
  message(chr, ': ', sum(ok), ' variants')
}
