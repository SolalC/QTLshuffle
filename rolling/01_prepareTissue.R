# =============================================================================
# 01  Tissue data and PEER per window
#
#   - donors of TISSUE with expression, covariates, a DIP and a genotype;
#   - expression of all genes (PEER is fitted on all of them), fixed covariates
#     (sex, platform, pcr, PC1-5 from the GTEx covariate file) and DIP;
#   - one window per centre in CENTRES_H (24 centres, 1 h apart): donors whose
#     DIP lies within WINDOW_WIDTH_H / 2 of the centre, wrapping around
#     midnight; centres with fewer than MIN_WINDOW_N donors are skipped;
#   - covariates of each window: fixed covariates plus peerK(n) PEER factors
#     recomputed on the window donors (15 if n < 150, as GTEx).
#
# Run:     sbatch 01_prepareTissue.sh
# Usage:   Rscript 01_prepareTissue.R <n cores>
# Outputs: Results_rolling/data/tissue.rds
#          Results_rolling/w<W>h/windowCovariates.rds: window table (centre, donor
#            rows, n, nPEER) and, per centre, covariates (donors x [fixed, PEER])
#          Results_rolling/w<W>h/windows.csv
# =============================================================================

.a <- commandArgs(FALSE); .f <- sub('^--file=', '', .a[grep('^--file=', .a)])
source(file.path(if (length(.f)) dirname(normalizePath(.f)) else '.', '00_config.R'))
library(parallel)
nCore <- as.numeric(commandArgs(trailingOnly = TRUE)[1])

tis <- loadTissue(TISSUE)
message(TISSUE, ': ', length(tis$donors), ' donors, ', ncol(tis$E), ' genes')
saveRDS(tis, TISSUE_FILE, compress = FALSE)

windows <- tibble(centreIdx = seq_along(CENTRES_H), centre = CENTRES_H) %>%
  mutate(rows = map(centre, \(c) windowRows(tis$radian, c)),
         n = lengths(rows), nPEER = peerK(n), used = n >= MIN_WINDOW_N)
print(select(windows, centre, n, nPEER, used), n = Inf)

# PEER per window, saved with the donor IDs as row names; re-used when the
# script is run again with the same window donors.
old <- if (file.exists(WINDOW_FILE)) readRDS(WINDOW_FILE) else NULL
if (!is.null(old) && identical(old$windows$rows, windows$rows)) {
  message('Re-using window covariates from ', WINDOW_FILE)
  cov <- old$cov
} else {
  cov <- mclapply(seq_len(nrow(windows)), \(i) {
    if (!windows$used[i]) return(NULL)
    windowCovariates(tis, windows$rows[[i]])
  }, mc.cores = nCore, mc.preschedule = FALSE)
  failed <- map_lgl(cov, \(x) inherits(x, 'try-error'))
  if (any(failed)) stop('PEER failed for centres ', paste(windows$centre[failed], collapse = ', '),
                        '\n', as.character(cov[[which(failed)[1]]]))
}

saveRDS(list(windows = windows, cov = cov), WINDOW_FILE)
fwrite(select(windows, centre, n, nPEER, used), paste0(WIDTH_DIR, 'windows.csv'))
