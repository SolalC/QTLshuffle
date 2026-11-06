# =============================================================================
# 05  24 h reference sample r: re-test of every window's independent signals
#
# One SLURM task per reference sample r = 1..N_REF. For each window centre c
# with nW donors (centres run in parallel over the task's cores):
#   1. nW donors drawn at random from all the tissue's donors, whatever their
#      DIP (window donors can be drawn; nShared records how many); the draw
#      depends only on (c, r, nW) through the seed;
#   2. covariates: fixed covariates plus peerK(nW) PEER factors recomputed on
#      those donors (same rule as the window);
#   3. every independent signal (gene, index variant) of centre c re-tested with
#      the same model; pairs with minor allele count < MIN_MAC are NA.
# The donors and covariates are saved and re-used if the task is run again: a
# reference sample does not depend on the signals, so a rerun of 04 (other
# thresholds, other clumping) needs no new PEER fits.
#
# Run:     sbatch 05_referenceRetest.sh   (array 1..N_REF, after 04)
# Usage:   Rscript 05_referenceRetest.R <reference sample r> <n cores>
# Outputs: Results/rolling/w<W>h/reference/ref<r>.csv.gz
#            one row per centre x signal: beta, se, p, nShared
#          Results/rolling/w<W>h/referenceCovariates/ref<r>.rds
#            list by centre: centre, donors, covariates (donors x [fixed, PEER])
# =============================================================================

.a <- commandArgs(FALSE); .f <- sub('^--file=', '', .a[grep('^--file=', .a)])
source(file.path(if (length(.f)) dirname(normalizePath(.f)) else '.', '00_config.R'))
library(parallel)
args <- commandArgs(trailingOnly = TRUE)
ref <- as.integer(args[1]); nCore <- as.integer(args[2])
if (ref > N_REF) stop('Reference sample ', ref, ' > N_REF')
outFile <- sprintf('%sref%04d.csv.gz', REF_DIR, ref)
covFile <- sprintf('%sref%04d.rds', REF_COV_DIR, ref)

win <- readRDS(WINDOW_FILE)$windows
tis <- readRDS(TISSUE_FILE)
used <- which(win$used)

# Reference donors and covariates (PEER) for every centre, computed once and saved.
old <- if (file.exists(covFile)) readRDS(covFile) else NULL
reuse <- !is.null(old) && identical(map_dbl(old, 'centre'), as.numeric(win$centre[used])) &&
  identical(map_int(old, \(x) length(x$donors)), as.integer(win$n[used]))
if (reuse) {
  refCov <- old
  message('Re-using reference covariates from ', covFile)
} else {
  t0 <- Sys.time()
  refCov <- mclapply(used, \(i) {
    rows <- referenceSample(length(tis$donors), win$n[i], i, ref)
    list(centre = win$centre[i], donors = tis$donors[rows], cov = windowCovariates(tis, rows))
  }, mc.cores = nCore, mc.preschedule = FALSE)
  failed <- map_lgl(refCov, \(x) inherits(x, 'try-error'))
  if (any(failed)) stop('PEER failed for centres ', paste(win$centre[used][failed], collapse = ', '),
                        '\n', as.character(refCov[[which(failed)[1]]]))
  saveRDS(refCov, covFile)
  message(length(used), ' PEER fits in ', round(as.numeric(difftime(Sys.time(), t0, units = 'mins')), 1), ' min')
}

signals <- fread(paste0(CALL_DIR, 'independent.csv.gz'))
if (!nrow(signals)) {
  message('No independent signals; nothing to re-test')
  fwrite(tibble(), outFile)
  quit(save = 'no')
}
geno <- openGenotypes(tis$donors)
G <- readGenotypes(geno, seq_along(tis$donors), match(unique(signals$variant), geno$bim$id))

res <- map_dfr(refCov, \(rc) {
  s <- filter(signals, centre == rc$centre)
  if (!nrow(s)) return(NULL)
  rows  <- match(rc$donors, tis$donors)
  rowsW <- win$rows[[match(rc$centre, win$centre)]]
  pairTest(tis$E[rows, s$gene, drop = FALSE], G[rows, s$variant, drop = FALSE], rc$cov) %>%
    mutate(centre = rc$centre, ref = ref, nShared = length(intersect(rows, rowsW)),
           gene = s$gene, variant = s$variant, .before = 1)
})
fwrite(res, outFile)
message(nrow(res), ' signal x centre re-tests')
