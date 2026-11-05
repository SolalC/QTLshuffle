# =============================================================================
# 05  24 h reference samples: re-test of each window's independent signals
#
# For window centre c with nW donors and reference sample r = 1..N_REF:
#   1. nW donors drawn at random from all the tissue's donors, whatever their
#      DIP (window donors can be drawn; nShared records how many);
#   2. covariates: fixed covariates plus peerK(nW) PEER factors recomputed on
#      those donors (same rule as the window);
#   3. every independent signal (gene, index variant) of centre c re-tested with
#      the same model; pairs with minor allele count < MIN_MAC in the sample
#      are NA.
# The donors and covariates (fixed + PEER) of every reference sample are saved
# per block and re-used if the task is run again: a reference sample depends
# only on the centre, nW and its seed, not on the signals, so a rerun of 04
# (other thresholds, other clumping) does not need new PEER fits.
#
# Run:     sbatch 05_referenceRetest.sh   (array 1..24 x N_REF / REF_BLOCK, after 04)
# Usage:   Rscript 05_referenceRetest.R <task> <n cores>
#          task = (centre index - 1) * (N_REF / REF_BLOCK) + block
# Outputs: Results_rolling/w<W>h/reference/c<centre>_block<b>.csv.gz
#          (one row per reference sample x signal)
#          Results_rolling/w<W>h/referenceCovariates/c<centre>_block<b>.rds: list,
#          one element per reference sample: ref, donors, covariates
#          (donors x [fixed covariates, PEER factors])
# =============================================================================

.a <- commandArgs(FALSE); .f <- sub('^--file=', '', .a[grep('^--file=', .a)])
source(file.path(if (length(.f)) dirname(normalizePath(.f)) else '.', '00_config.R'))
library(parallel)
args <- commandArgs(trailingOnly = TRUE)
task <- as.integer(args[1]); nCore <- as.integer(args[2])
nBlock    <- ceiling(N_REF / REF_BLOCK)
centreIdx <- (task - 1) %/% nBlock + 1
block     <- (task - 1) %% nBlock + 1
centre    <- CENTRES_H[centreIdx]
refs      <- ((block - 1) * REF_BLOCK + 1):min(block * REF_BLOCK, N_REF)
outFile   <- sprintf('%sc%02d_block%d.csv.gz', REF_DIR, centre, block)
covFile   <- sprintf('%sc%02d_block%d.rds', REF_COV_DIR, centre, block)

win <- readRDS(WINDOW_FILE)$windows
if (!win$used[centreIdx]) {
  message('Centre ', centre, ' h has fewer than ', MIN_WINDOW_N, ' donors: skipped')
  quit(save = 'no')
}
tis <- readRDS(TISSUE_FILE)
rowsW <- win$rows[[centreIdx]]
nW <- length(rowsW)

# Reference samples and their covariates (PEER), computed once and saved.
refCov <- if (file.exists(covFile)) readRDS(covFile) else NULL
if (is.null(refCov) || !identical(map_int(refCov, 'ref'), refs) ||
    !all(map_int(refCov, \(x) length(x$donors)) == nW)) {
  t0 <- Sys.time()
  refCov <- mclapply(refs, \(ref) {
    rows <- referenceSample(length(tis$donors), nW, centreIdx, ref)
    list(ref = ref, donors = tis$donors[rows], cov = windowCovariates(tis, rows))
  }, mc.cores = nCore, mc.preschedule = FALSE)
  failed <- map_lgl(refCov, \(x) inherits(x, 'try-error'))
  if (any(failed)) stop('PEER failed for references ', paste(refs[failed], collapse = ', '), '\n',
                        as.character(refCov[[which(failed)[1]]]))
  saveRDS(refCov, covFile)
  message(length(refs), ' reference PEER fits in ',
          round(as.numeric(difftime(Sys.time(), t0, units = 'mins')), 1), ' min')
} else message('Re-using reference covariates from ', covFile)

signals <- fread(paste0(CALL_DIR, 'independent.csv.gz')) %>% filter(centre == !!centre)
if (!nrow(signals)) {
  message('Centre ', centre, ' h: no independent signals; nothing to re-test')
  fwrite(tibble(), outFile)
  quit(save = 'no')
}
G <- dosage(readRDS(paste0(CALL_DIR, 'indexGenotypes.rds'))[, signals$variant, drop = FALSE])
stopifnot(identical(rownames(G), tis$donors))
Y <- tis$E[, signals$gene, drop = FALSE]
message('Centre ', centre, ' h: ', nrow(signals), ' signals, nW = ', nW, ', references ',
        min(refs), '-', max(refs))

res <- map_dfr(refCov, \(rc) {
  rows <- match(rc$donors, tis$donors)
  pairTest(Y[rows, , drop = FALSE], G[rows, , drop = FALSE], covariateBasis(rc$cov)) %>%
    mutate(centre = centre, ref = rc$ref, nShared = length(intersect(rows, rowsW)),
           gene = signals$gene, variant = signals$variant, .before = 1)
})
fwrite(res, outFile)
