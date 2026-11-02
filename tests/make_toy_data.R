# =============================================================================
# Toy GTEx-like inputs for tests/run_toy_test.sh. Nothing here is real data.
#
# Two tissues (ToyA: 560 donors, ToyB: 520 donors), 2 chromosomes x 300
# variants, 2 x 150 genes. DIP is uniform on the circle. Gene architectures:
#   genes  1-10  phase QTL     (no 24 h mean effect; 3 h shift per allele)
#   genes 11-20  amplitude QTL (no 24 h mean effect)
#   genes 21-30  mean eQTL     (time-invariant)
#   genes 31-40  rhythmic, no genetic effect
#   genes 41-300 null (enough genes for 10 latent factors)
# Phase and amplitude QTL should give vanishing hits; mean eQTL should not.
#
# Usage: Rscript make_toy_data.R <output dir>
# =============================================================================

suppressPackageStartupMessages({ library(tidyverse); library(data.table) })
set.seed(20261102)
out <- commandArgs(trailingOnly = TRUE)[1]
for (d in c('expr', 'cov', 'harmonic', 'chiral'))
  dir.create(file.path(out, d), recursive = TRUE, showWarnings = FALSE)

nDonor <- 600; nVarChr <- 300; nGeneChr <- 150
donors <- sprintf('GTEX-T%03d', seq_len(nDonor))
phi <- runif(nDonor, -pi, pi); names(phi) <- sub('GTEX-', '', donors)
save(phi, file = file.path(out, 'DIP.RData'))

# Genotypes (genio layout: X = variants x WGS samples, bim):
bim <- map_dfr(1:2, \(c) tibble(chr = as.character(c), id = sprintf('rs%d%04d', c, 1:nVarChr),
                                pos = sort(sample(1:3e6, nVarChr))))
maf <- runif(nrow(bim), 0.15, 0.5)
X <- t(sapply(maf, \(p) rbinom(nDonor, 2, p)))
X[sample(length(X), 200)] <- NA
dimnames(X) <- list(bim$id, paste0(donors, '-0003-SM-TOY'))
saveRDS(list(X = X, bim = as.data.frame(bim)), file.path(out, 'genotypes.rds'))

genes <- map_dfr(1:2, \(c) tibble(`#chr` = paste0('chr', c),
  start = round(seq(2e5, 2.8e6, length.out = nGeneChr)), gene_id = sprintf('ENSG%02d%03d.1', c, 1:nGeneChr))) %>%
  mutate(end = start + 1, k = rep(seq_len(nGeneChr), 2), .after = start)

simulateExpr <- function(donorIdx) {
  t <- phi[donorIdx]
  map(seq_len(nrow(genes)), \(k) {
    chr <- as.integer(sub('chr', '', genes$`#chr`[k]))
    onChr <- which(bim$chr == chr)
    causal <- onChr[which.min(abs(bim$pos[onChr] - genes$start[k]))]
    g <- X[causal, donorIdx]; g[is.na(g)] <- 1
    peak <- runif(1, -pi, pi)
    y <- case_when(
      k %in% 1:10  ~ 2.0 * cos(t - peak - (pi / 4) * g),
      k %in% 11:20 ~ 1.2 * g * cos(t - peak),
      k %in% 21:30 ~ 1.0 * g,
      k %in% 31:40 ~ 0.8 * cos(t - peak),
      TRUE ~ 0)
    y + rnorm(length(t))
  })
}

writeTissue <- function(name, n) {
  idx <- sort(sample(nDonor, n))
  E <- simulateExpr(idx)
  expr <- bind_cols(genes %>% select(`#chr`, start, end, gene_id),
                    as_tibble(do.call(rbind, E), .name_repair = \(x) donors[idx]))
  fwrite(expr, file.path(out, 'expr', paste0(name, '.v10.normalized_expression.bed.gz')), sep = '\t')
  covRows <- c('sex', 'platform', 'pcr', paste0('PC', 1:5), paste0('InferredCov', 1:10))
  cov <- matrix(rnorm(length(covRows) * n), nrow = length(covRows))
  cov[1, ] <- sample(1:2, n, TRUE); cov[2:3, ] <- sample(0:1, 2 * n, TRUE)
  fwrite(cbind(data.table(ID = covRows), as.data.table(`colnames<-`(cov, donors[idx]))),
         file.path(out, 'cov', paste0(name, '.v10.covariates.txt')), sep = '\t')
  fwrite(tibble(gene = genes$gene_id,
                p.value = ifelse(genes$k %in% c(1:20, 31:40), 1e-12, 0.5)),
         file.path(out, 'harmonic', paste0(name, '_LRT_harmonic_linear.csv')))
}
writeTissue('ToyA', 560)
writeTissue('ToyB', 520)
# One gene stands in for a CHIRAL input gene and should be dropped by 02:
writeLines('ENSG02150', file.path(out, 'chiral', 'CHIRAL_genes.txt'))
