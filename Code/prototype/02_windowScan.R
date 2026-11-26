# =============================================================================
# 02  Window scan with the permuted-DIP null (one block of DIP permutations per
#     SLURM array task)
#
# For each DIP permutation (0 = observed; 1..N_PERM = DIP shuffled across
# donors) and each of the N_WINDOWS 3 h DIP windows with >= MIN_WINDOW_N donors:
#   1. PEER recomputed on the window donors (all genes);
#   2. each selected gene scanned as tensorQTL map_cis does (NPERM_GENO
#      permutations, beta-approximated p); eGenes: Storey q < GENE_FDR;
#   3. 24 h reference: as many donors as in the window, drawn from all donors
#      (about 1/8 overlap); PEER recomputed there; lead variants re-tested;
#   4. vanishing hit: reference p > VANISH_P and z test p < HET_P.
#
# Run:     sbatch 02_windowScan.sh          (array 1..ceiling((N_PERM+1)/PERM_BLOCK))
# Usage:   Rscript 02_windowScan.R <block> <n cores>
# Outputs: Results/prototype/windowCounts/block<b>.csv  (one row per perm x window)
#          Results/prototype/windowHits/observed.csv    (block 1 only: perm 0)
# =============================================================================

.a <- commandArgs(FALSE); .f <- sub('^--file=', '', .a[grep('^--file=', .a)])
source(file.path(if (length(.f)) dirname(normalizePath(.f)) else '.', '00_config.R'))
library(parallel)

args  <- commandArgs(trailingOnly = TRUE)
block <- as.numeric(args[1])
nCore <- as.numeric(args[2])
perms <- seq((block - 1) * PERM_BLOCK, min(block * PERM_BLOCK - 1, N_PERM))
if (!length(perms) || min(perms) > N_PERM) stop('Block ', block, ' is beyond N_PERM')
message('Permutations ', min(perms), '-', max(perms))
d <- readRDS(DATA_FILE)

runPermutation <- function(perm) {
  t0  <- Sys.time()
  win <- dipWindow(permutedRadian(d$tis$radian, perm))
  out <- map(seq_len(N_WINDOWS), \(w) {
    rowsW <- which(win == w)
    if (length(rowsW) < MIN_WINDOW_N) return(list(counts = tibble(nW = length(rowsW)), hits = NULL))
    rowsR <- referenceSample(length(win), length(rowsW), perm, w)
    scanWindow(d, rowsW, rowsR, seedBase = SEED + 1e5 * w + 100 * perm)
  })
  counts <- map_dfr(seq_len(N_WINDOWS), \(w) mutate(out[[w]]$counts, window = w)) %>%
    mutate(tissue = TISSUE, perm = perm, .before = 1)
  hits <- if (perm == 0)
    map_dfr(seq_len(N_WINDOWS), \(w) if (!is.null(out[[w]]$hits)) mutate(out[[w]]$hits, window = w))
  message('perm ', perm, ': ', round(as.numeric(difftime(Sys.time(), t0, units = 'mins')), 2), ' min')
  list(counts = counts, hits = hits)
}

res <- mclapply(perms, runPermutation, mc.cores = nCore, mc.preschedule = FALSE)
failed <- map_lgl(res, \(r) inherits(r, 'try-error'))
if (any(failed)) stop('Permutations failed: ', paste(perms[failed], collapse = ', '),
                      '\n', as.character(res[[which(failed)[1]]]))

fwrite(map_dfr(res, 'counts'), paste0(COUNT_DIR, 'block', block, '.csv'))
if (0 %in% perms) fwrite(res[[which(perms == 0)]]$hits, paste0(HIT_DIR, 'observed.csv'))
