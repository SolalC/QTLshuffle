# =============================================================================
# 08  Selection-bias control, pass B: re-tests and calibration, one centre
#
# For window centre c (one SLURM task per centre):
#   1. pairs: the window's independent signals (04) and the pseudo-window
#      signals of references 1..N_NULL (07);
#   2. beta of every pair in each reference k = 1..N_NULL of centre c (its
#      donors and saved covariates), with MatrixEQTL (pairBeta, the model of
#      the scan); the beta of a pseudo-window signal in its own reference is
#      recomputed as a check (must equal the 07 beta);
#   3. zRef = (beta - mean) / SD of the comparison betas, and pEmp (two-sided,
#      as 06): window signals against all N_NULL references, pseudo-window
#      signals of reference k against the other N_NULL - 1;
#   4. calibrated p of each window signal: within its |z| bin (Z_BINS, z =
#      beta / se in the window),
#        pCal = (1 + #{pseudo-window signals with |zRef| >= its |zRef|}) / (1 + #{pseudo-window signals})
#      i.e. how often a signal selected the same way, of the same strength,
#      but without time structure deviates this much from its references;
#      qCal = BH of pCal over the window's signals; time dependent: qCal < TIME_FDR;
#   5. window level: % of window signals with pEmp < 0.05 against the same %
#      in each pseudo-window (the % produced by selection alone);
#      pWindow = (1 + #{pseudo-windows with % >= window %}) / (1 + N_NULL).
#
# Run:     sbatch 08_nullCalibrate.sh   (array 1..24 centres, after every 07 task)
# Usage:   Rscript 08_nullCalibrate.R <centre index> <n cores>
# Outputs: Results/rolling/w<W>h/null/calibration/
#            c<cc>_window.csv.gz   window signals: zW, zBin, nCmp, refMeanK, refSDK,
#                                  zRefK, pEmpK, pCal, qCal, timeDependent
#            c<cc>_null.csv.gz     pseudo-window signals: ref, zW, zBin, nCmp, zRef, pEmp, betaSelf
#            c<cc>_nullByRef.csv   % pEmp < 0.05 per pseudo-window
#            c<cc>_summary.csv     window %, null % (mean, SD, range), pWindow, time-dependent counts
#            c<cc>_betas.rds       pairs x N_NULL reference betas (re-used if present)
# =============================================================================

.a <- commandArgs(FALSE); .f <- sub('^--file=', '', .a[grep('^--file=', .a)])
source(file.path(if (length(.f)) dirname(normalizePath(.f)) else '.', '00_config.R'))
args <- commandArgs(trailingOnly = TRUE)
centreIdx <- as.integer(args[1]); nCore <- as.integer(args[2])
centre <- CENTRES_H[centreIdx]
out <- \(what) sprintf('%sc%02d_%s', NULL_CAL_DIR, centre, what)

win <- readRDS(WINDOW_FILE)$windows
if (!win$used[centreIdx]) { message('Centre ', centre, ' h skipped (< ', MIN_WINDOW_N, ' donors)'); quit(save = 'no') }

# --- 1. Pairs -------------------------------------------------------------------------
cc <- centre   # (data.table: avoid the column name)
window <- fread(paste0(CALL_DIR, 'independent.csv.gz'))[centre == cc]
nullFiles <- sprintf('%sref%04d.csv.gz', NULL_SCAN_DIR, seq_len(N_NULL))
if (!all(file.exists(nullFiles))) stop(sum(!file.exists(nullFiles)), ' 07 outputs missing, e.g. ', nullFiles[!file.exists(nullFiles)][1])
pseudo <- rbindlist(map(nullFiles, \(f) { x <- fread(f); if (nrow(x)) x[centre == cc] else NULL }))
# Nothing to calibrate: no window signal, or no pseudo-window signal at this centre.
# Only the summary is written (window file absent; 09 skips the centre).
if (!nrow(window) || !nrow(pseudo)) {
  message('Centre ', centre, ' h: ', nrow(window), ' window and ', nrow(pseudo),
          ' pseudo-window signals: nothing to calibrate')
  fwrite(data.table(centre = centre, nW = win$n[centreIdx], windowSignals = nrow(window),
                    pseudoWindows = uniqueN(pseudo$ref), pseudoSignalsMean = nrow(pseudo) / N_NULL,
                    pctWindowPEmp05 = NA_real_, pctNullMean = NA_real_, pctNullSD = NA_real_,
                    pctNullMin = NA_real_, pctNullMax = NA_real_, pWindow = NA_real_,
                    timeDependentCal = 0L, untestedCal = nrow(window)), out('summary.csv'))
  quit(save = 'no')
}
window[, gene := as.character(gene)]; pseudo[, gene := as.character(gene)]
pairs <- unique(rbind(window[, .(gene, variant, chr)], pseudo[, .(gene, variant, chr)]))
pairs[, id := .I]
message('Centre ', centre, ' h: ', nrow(window), ' window signals, ', nrow(pseudo), ' pseudo-window signals in ',
        uniqueN(pseudo$ref), ' references, ', nrow(pairs), ' distinct pairs')

# --- 2. Betas of every pair in references 1..N_NULL -------------------------------------
betaFile <- out('betas.rds')
old <- if (file.exists(betaFile)) readRDS(betaFile) else NULL
if (!is.null(old) && identical(old$pairs[, .(gene, variant)], pairs[, .(gene, variant)]) && ncol(old$B) == N_NULL) {
  B <- old$B
  message('Re-using reference betas from ', betaFile)
} else {
  tis  <- readRDS(TISSUE_FILE)
  geno <- openGenotypes(tis$donors)
  t0 <- Sys.time()
  cols <- parallel::mclapply(seq_len(N_NULL), \(k) {
    rc <- keep(readRDS(sprintf('%sref%04d.rds', REF_COV_DIR, k)), \(x) x$centre == centre)[[1]]
    rows <- match(rc$donors, tis$donors)
    b <- rep(NA_real_, nrow(pairs))
    for (chrom in unique(pairs$chr)) {
      p <- pairs[chr == chrom]
      r <- pairBeta(tis, geno, rows, rc$cov, p, chrom)
      b[p$id] <- r$beta
    }
    b
  }, mc.cores = nCore, mc.preschedule = FALSE)
  failed <- map_lgl(cols, \(x) inherits(x, 'try-error'))
  if (any(failed)) stop('Re-test failed for ', sum(failed), ' references: ', as.character(cols[[which(failed)[1]]]))
  B <- do.call(cbind, cols)
  saveRDS(list(pairs = pairs, B = B), betaFile)
  message(N_NULL, ' references re-tested in ', round(as.numeric(difftime(Sys.time(), t0, units = 'mins')), 1), ' min')
}

# --- 3. zRef and pEmp ---------------------------------------------------------------------
# b: beta per signal; M: signals x comparison references (rows aligned)
compare <- function(b, M) {
  n <- rowSums(!is.na(M))
  m <- rowMeans(M, na.rm = TRUE)
  s <- sqrt(rowSums((M - m)^2, na.rm = TRUE) / (n - 1))
  data.table(nCmp = n, refMean = m, refSD = s, zRef = (b - m) / s,
             pEmp = (1 + rowSums(abs(M - m) >= abs(b - m), na.rm = TRUE)) / (1 + n))
}
window[pairs, id := i.id, on = c('gene', 'variant')]
pseudo[pairs, id := i.id, on = c('gene', 'variant')]

w <- compare(window$beta, B[window$id, , drop = FALSE])
window <- cbind(window, setnames(w, c('nCmp', 'refMeanK', 'refSDK', 'zRefK', 'pEmpK')))
pseudo <- rbindlist(map(split(pseudo, by = 'ref'), \(x)        # reference k: compared with the other N_NULL - 1
  cbind(x, compare(x$beta, B[x$id, -x$ref[1], drop = FALSE]))))
pseudo[, betaSelf := B[cbind(id, ref)]]
selfDiff <- max(abs(pseudo$betaSelf - pseudo$beta), na.rm = TRUE)
message('Check: pseudo-window beta recomputed in its own reference, max |difference| = ', signif(selfDiff, 3))
if (selfDiff > 1e-6) warning('pairBeta does not reproduce the 07 scan beta (max difference ', signif(selfDiff, 3), ')')

# --- 4. Calibrated p within |z| bins -----------------------------------------------------
window[, `:=`(zW = beta / se, zBin = cut(abs(beta / se), Z_BINS))]
pseudo[, `:=`(zW = beta / se, zBin = cut(abs(beta / se), Z_BINS))]
nullAbs <- split(abs(pseudo$zRef[!is.na(pseudo$zRef)]), pseudo$zBin[!is.na(pseudo$zRef)])
window[, nNullBin := map_int(as.character(zBin), \(b) length(nullAbs[[b]]))]
window[, pCal := map2_dbl(as.character(zBin), abs(zRefK), \(b, z) {
  x <- nullAbs[[b]]
  if (is.na(z) || !length(x)) NA_real_ else (1 + sum(x >= z)) / (1 + length(x))
})]
window[, qCal := p.adjust(pCal, 'BH')]
window[, timeDependent := qCal < TIME_FDR]

# --- 5. Window level ------------------------------------------------------------------------
nullByRef <- pseudo[, .(signals = .N, pctPEmp05 = 100 * mean(pEmp < 0.05, na.rm = TRUE)), by = ref][order(ref)]
pctWindow <- 100 * mean(window$pEmpK < 0.05, na.rm = TRUE)
summary <- data.table(
  centre = centre, nW = win$n[centreIdx], windowSignals = nrow(window), pseudoWindows = nrow(nullByRef),
  pseudoSignalsMean = mean(nullByRef$signals),
  pctWindowPEmp05 = pctWindow, pctNullMean = mean(nullByRef$pctPEmp05), pctNullSD = sd(nullByRef$pctPEmp05),
  pctNullMin = min(nullByRef$pctPEmp05), pctNullMax = max(nullByRef$pctPEmp05),
  pWindow = (1 + sum(nullByRef$pctPEmp05 >= pctWindow)) / (1 + nrow(nullByRef)),
  timeDependentCal = sum(window$timeDependent, na.rm = TRUE), untestedCal = sum(is.na(window$pCal)))
print(summary)
print(window[, .(signals = .N, pctPEmp05 = round(100 * mean(pEmpK < 0.05, na.rm = TRUE), 1),
                 pctTimeDependentCal = round(100 * mean(timeDependent, na.rm = TRUE), 1), nullSignals = nNullBin[1]),
             by = zBin][order(zBin)])

fwrite(window[, .(centre, gene, chr, variant, signal, nW, beta, se, pNominal, qval, zW, zBin, nNullBin,
                  nCmp, refMeanK, refSDK, zRefK, pEmpK, pCal, qCal, timeDependent)], out('window.csv.gz'))
fwrite(pseudo[, .(centre, ref, gene, chr, variant, signal, beta, se, pNominal, qval, zW, zBin,
                  nCmp, refMean, refSD, zRef, pEmp, betaSelf)], out('null.csv.gz'))
fwrite(nullByRef, out('nullByRef.csv'))
fwrite(summary, out('summary.csv'))
