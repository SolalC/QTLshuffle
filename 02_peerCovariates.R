# =============================================================================
# 02  Covariates for every window and 24 h reference sample (one tissue x
#     permutation block per SLURM array task)
#
# For each permutation (0 = observed DIP; 1..N_PERM = DIP shuffled across
# donors) and each of the N_WINDOWS 3 h DIP windows:
#   - window donors: DIP in the window (skipped if fewer than MIN_WINDOW_N);
#   - 24 h reference: as many donors as in the window, drawn at random from all
#     the tissue's donors whatever their DIP (about 1/8 are window donors);
#   - covariates of each set: sex, platform, pcr, PC1-5 and N_PEER PEER factors
#     recomputed on that set's expression.
#
# Run:     sbatch 02_peerCovariates.sh      (array: N_TISSUES x blocks; after 01)
# Usage:   Rscript 02_peerCovariates.R <task index> <n cores>
#          task -> tissue = (task-1) %/% nBlock + 1, block = (task-1) %% nBlock + 1,
#          nBlock = ceiling((N_PERM + 1) / PERM_BLOCK)
# Outputs: Results/peerCovariates/<tissue>_block<b>.csv.gz, long format:
#          perm, window, set (window / reference), donor, one column per covariate
# =============================================================================

.a <- commandArgs(FALSE); .f <- sub('^--file=', '', .a[grep('^--file=', .a)])
source(file.path(if (length(.f)) dirname(normalizePath(.f)) else '.', '00_config.R'))
library(parallel)

args  <- commandArgs(trailingOnly = TRUE)
task  <- as.numeric(args[1])
nCore <- as.numeric(args[2])

tissues   <- fread(paste0(TISSUE_DIR, 'selectedTissues.csv'))
nBlock    <- ceiling((N_PERM + 1) / PERM_BLOCK)
tissueIdx <- (task - 1) %/% nBlock + 1
block     <- (task - 1) %% nBlock + 1
if (tissueIdx > nrow(tissues)) stop('Task ', task, ' is beyond the last tissue')
perms <- seq((block - 1) * PERM_BLOCK, min(block * PERM_BLOCK - 1, N_PERM))
tissueName <- tissues$tissue[tissueIdx]
message('Tissue: ', tissueName, '; permutations ', min(perms), '-', max(perms))
tis <- readRDS(paste0(TISSUE_DIR, tissueName, '.rds'))

covTable <- function(rows, perm, w, set) {
  C <- subsetCovariates(tis, rows)
  data.table(perm = perm, window = w, set = set, donor = tis$donors[rows], C)
}

runPermutation <- function(perm) {
  t0  <- Sys.time()
  win <- dipWindow(permutedRadian(tis$radian, perm, tissueIdx))
  out <- map(seq_len(N_WINDOWS), \(w) {
    rowsW <- which(win == w)
    if (length(rowsW) < MIN_WINDOW_N) return(NULL)
    rowsR <- referenceSample(length(win), length(rowsW), perm, tissueIdx, w)
    rbind(covTable(rowsW, perm, w, 'window'), covTable(rowsR, perm, w, 'reference'),
          fill = TRUE)
  })
  message('perm ', perm, ': ', round(as.numeric(difftime(Sys.time(), t0, units = 'mins')), 2), ' min')
  rbindlist(out, fill = TRUE)
}

res <- mclapply(perms, runPermutation, mc.cores = nCore, mc.preschedule = FALSE)
failed <- map_lgl(res, \(r) inherits(r, 'try-error'))
if (any(failed)) stop('Permutations failed: ', paste(perms[failed], collapse = ', '),
                      '\n', as.character(res[[which(failed)[1]]]))

# Covariates dropped in a subset (constant there) are left empty (NA):
fwrite(rbindlist(res, fill = TRUE), paste0(PEER_DIR, tissueName, '_block', block, '.csv.gz'))
