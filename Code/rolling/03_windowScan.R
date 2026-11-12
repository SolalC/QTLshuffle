# =============================================================================
# 03  cis-eQTL scan of one window centre x one chromosome (MatrixEQTL)
#
# Window donors only. For every gene on the chromosome (TSS = `end` of the
# expression .bed), every cis variant within CIS_WINDOW of the TSS with minor
# allele count >= MIN_MAC among the window donors (missing genotypes
# mean-imputed) is tested with MatrixEQTL's linear model:
#   y ~ g + sex + platform + pcr + PC1-5 + peerK(n) PEER factors
# (window covariates from 01). Nominal p-values only: no permutations. Pairs
# with p < PAIR_P_STORE are stored, with the number of cis tests, for the
# FDR in 04. Genotypes are read from the plink file of 02a (BEDMatrix).
#
# Run:     sbatch 03_windowScan.sh       (array 1..24 x length(CHROMOSOMES))
# Usage:   Rscript 03_windowScan.R <task>
#          task = (centre index - 1) * length(CHROMOSOMES) + chromosome index
# Outputs: Results/rolling/w<W>h/scan/c<centre>_chr<c>_{pairs,genes}.csv.gz
#            pairs: gene, variant, beta, se, t, pNominal (p < PAIR_P_STORE)
#            genes: gene, nVar (cis variants tested), minP; nTests in every row
# =============================================================================

.a <- commandArgs(FALSE); .f <- sub('^--file=', '', .a[grep('^--file=', .a)])
source(file.path(if (length(.f)) dirname(normalizePath(.f)) else '.', '00_config.R'))
task <- as.integer(commandArgs(trailingOnly = TRUE)[1])
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
rows <- win$windows$rows[[centreIdx]]
geno <- openGenotypes(tis$donors)

t0 <- Sys.time()
sc <- scanChromosome(tis, geno, rows, win$cov[[centreIdx]], chrom)   # 00_functions.R
message('Centre ', centre, ' h, chr', chrom, ': ', length(rows), ' donors, ', sc$nCov, ' covariates, ',
        sc$nVariants, ' variants with MAC >= ', MIN_MAC)
message('MatrixEQTL: ', sc$nTests, ' cis tests in ',
        round(as.numeric(difftime(Sys.time(), t0, units = 'mins')), 1), ' min')

fwrite(mutate(sc$genes, centre = centre, chr = chrom, .before = 1) %>% mutate(nTests = sc$nTests),
       outFile('genes'))
fwrite(mutate(sc$pairs, centre = centre, .before = 1), outFile('pairs'))
