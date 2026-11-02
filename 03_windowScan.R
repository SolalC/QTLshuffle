# =============================================================================
# 03  Window scan with permuted-DIP null (one tissue x permutation block per task)
#
# For each permutation (0 = observed DIP; 1..N_PERM = DIP shuffled across donors)
# and each of the N_WINDOWS 3 h DIP windows:
#   1. PEER recomputed on the window donors; cis-eQTL scan (all genes, all cis
#      variants with MAC >= MIN_MAC in the window);
#   2. eGenes: lead p x Meff, BH across genes, q < GENE_FDR;
#   3. PEER recomputed on the out-of-window donors; each eGene's lead variant
#      re-tested there;
#   4. vanishing hit: out-of-window p > VANISH_P and window vs out-of-window
#      difference z test p < HET_P.
#
# Run:     sbatch 03_windowScan.sh          (array: N_TISSUES x blocks; after 02)
# Usage:   Rscript 03_windowScan.R <task index> <n cores>
#          task -> tissue = (task-1) %/% nBlock + 1, block = (task-1) %% nBlock + 1,
#          nBlock = ceiling((N_PERM + 1) / PERM_BLOCK)
# Outputs: Results/windowCounts/<tissue>_block<b>.csv  (one row per perm x window)
#          Results/windowHits/<tissue>_observed.csv    (block 1 only: hits, perm 0)
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

tis  <- readRDS(paste0(TISSUE_DIR, tissueName, '.rds'))
geno <- map(set_names(names(tis$cis)), \(chr) list(G = readRDS(paste0(GENO_DIR, chr, '.rds'))$G))

runPermutation <- function(perm) {
  t0  <- Sys.time()
  win <- dipWindow(permutedRadian(tis$radian, perm, tissueIdx))
  out <- map(seq_len(N_WINDOWS), \(w) scanWindow(tis, geno, win == w))
  counts <- map_dfr(seq_len(N_WINDOWS), \(w) out[[w]]$counts %>% mutate(window = w)) %>%
    mutate(tissue = tissueName, perm = perm, .before = 1)
  hits <- if (perm == 0)
    map_dfr(seq_len(N_WINDOWS), \(w) if (!is.null(out[[w]]$hits)) mutate(out[[w]]$hits, window = w)) %>%
      mutate(tissue = tissueName, .before = 1)
  message('perm ', perm, ': ', round(as.numeric(difftime(Sys.time(), t0, units = 'mins')), 2), ' min')
  list(counts = counts, hits = hits)
}

res <- mclapply(perms, runPermutation, mc.cores = nCore, mc.preschedule = FALSE)
failed <- map_lgl(res, \(r) inherits(r, 'try-error'))
if (any(failed)) stop('Permutations failed: ', paste(perms[failed], collapse = ', '),
                      '\n', as.character(res[[which(failed)[1]]]))

fwrite(map_dfr(res, 'counts'), paste0(COUNT_DIR, tissueName, '_block', block, '.csv'))
if (0 %in% perms) fwrite(res[[which(perms == 0)]]$hits, paste0(HIT_DIR, tissueName, '_observed.csv'))
