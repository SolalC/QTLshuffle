# =============================================================================
# 03  cis-eQTL scan of one window centre x one chromosome
#
# For every gene on the chromosome (TSS from the expression .bed), cis variants
# within CIS_WINDOW of the TSS with minor allele count >= MIN_MAC among the
# window donors, tested with the window covariates from 01:
#   y ~ g + sex + platform + pcr + PC1-5 + peerK(n) PEER factors.
# Genotypes are read from the plink file of 02a (BEDMatrix).
# Gene level: lead variant, NPERM_GENO expression permutations (Freedman-Lane),
# beta-approximated p (tensorQTL). Pair level: every pair with nominal
# p < PAIR_P_STORE is kept for 04.
#
# Run:     sbatch 03_windowScan.sh       (array 1..24 x length(CHROMOSOMES))
# Usage:   Rscript 03_windowScan.R <task> <n cores>
#          task = (centre index - 1) * length(CHROMOSOMES) + chromosome index
# Outputs: Results/rolling/w<W>h/scan/c<centre>_chr<c>_{genes,pairs}.csv.gz
# =============================================================================

.a <- commandArgs(FALSE); .f <- sub('^--file=', '', .a[grep('^--file=', .a)])
source(file.path(if (length(.f)) dirname(normalizePath(.f)) else '.', '00_config.R'))
library(parallel)
args <- commandArgs(trailingOnly = TRUE)
task <- as.integer(args[1]); nCore <- as.integer(args[2])
centreIdx <- (task - 1) %/% length(CHROMOSOMES) + 1
chrom     <- CHROMOSOMES[(task - 1) %% length(CHROMOSOMES) + 1]
centre    <- CENTRES_H[centreIdx]
outFile   <- \(what) sprintf('%sc%02d_chr%d_%s.csv.gz', SCAN_DIR, centre, chrom, what)

win <- readRDS(WINDOW_FILE)
if (!win$windows$used[centreIdx]) {
  message('Centre ', centre, ' h has ', win$windows$n[centreIdx], ' donors (< ', MIN_WINDOW_N, '): skipped')
  fwrite(tibble(), outFile('genes')); fwrite(tibble(), outFile('pairs'))
  quit(save = 'no')
}
tis  <- readRDS(TISSUE_FILE)
geno <- openGenotypes(tis$donors)
rows <- win$windows$rows[[centreIdx]]
Q <- covariateBasis(win$cov[[centreIdx]])
message('Centre ', centre, ' h, chr', chrom, ': ', length(rows), ' donors, ',
        ncol(win$cov[[centreIdx]]), ' covariates')

genes <- tis$pos[sub('^chr', '', chr) == as.character(chrom)]
genes[, geneIdx := match(gene_id, colnames(tis$E))]
cols <- which(geno$bim$chr == as.character(chrom))
pos  <- geno$bim$pos[cols]
stopifnot(!is.unsorted(pos))
Gw <- readGenotypes(geno, rows, cols)          # window donors x chromosome variants
lo <- findInterval(genes$start - CIS_WINDOW - 1, pos) + 1
hi <- findInterval(genes$start + CIS_WINDOW, pos)
message(nrow(genes), ' genes, ', ncol(Gw), ' variants')

t0 <- Sys.time()
res <- mclapply(seq_len(nrow(genes)), \(i) {
  if (hi[i] < lo[i]) return(NULL)
  out <- geneScan(tis$E[rows, genes$geneIdx[i]], Gw[, lo[i]:hi[i], drop = FALSE], Q,
                  seed = SEED + 1e5 * centreIdx + genes$geneIdx[i])
  if (!is.null(out)) map(out, \(x) mutate(x, gene = genes$gene_id[i], .before = 1))
}, mc.cores = nCore, mc.preschedule = TRUE)
failed <- map_lgl(res, \(x) inherits(x, 'try-error'))
if (any(failed)) stop(sum(failed), ' genes failed: ', as.character(res[[which(failed)[1]]]))
res <- compact(res)
message(length(res), ' genes tested in ', round(as.numeric(difftime(Sys.time(), t0, units = 'mins')), 1), ' min')

fwrite(map_dfr(res, 'lead') %>% mutate(centre = centre, chr = chrom, .before = 1), outFile('genes'))
fwrite(map_dfr(res, 'pairs') %>% mutate(centre = centre, .before = 1), outFile('pairs'))
