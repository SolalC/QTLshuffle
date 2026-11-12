# =============================================================================
# 04  Significant pairs, eGenes and independent signals per window centre
#
# For each centre:
#   1. Benjamini-Hochberg q-value of every stored pair over all cis tests of
#      the window (the sum of MatrixEQTL's test counts over chromosomes, so
#      the same FDR as one MatrixEQTL run on the whole genome); significant
#      pair: q < PAIR_FDR; eGene: a gene with at least one significant pair.
#      Nominal p-values without permutations, so this is a pair-level FDR,
#      not the GTEx gene-level permutation threshold;
#   2. independent signals: each eGene's significant variants clumped with
#      plink2 --clump (r2 > CLUMP_R2 within CLUMP_KB kb; p1 = p2 = 1, i.e.
#      every significant variant can be an index variant). LD is computed by
#      plink2 on all the tissue's donors (the per-chromosome plink files of
#      02a), not the window donors: with ~90 donors, r2 between unlinked
#      variants averages ~1/90, above the 0.01 cut-off.
#
# Run:     sbatch 04_windowCall.sh      (after every 03 task)
# Usage:   Rscript 04_windowCall.R <n cores>
# Outputs: Results/rolling/w<W>h/calls/
#            genes.csv.gz            every gene tested, per centre (nVar, minP, eGene)
#            significantPairs.csv.gz pairs with q < PAIR_FDR
#            independent.csv.gz      plink2 index variants (one row per signal)
# =============================================================================

.a <- commandArgs(FALSE); .f <- sub('^--file=', '', .a[grep('^--file=', .a)])
source(file.path(if (length(.f)) dirname(normalizePath(.f)) else '.', '00_config.R'))
nCore <- as.integer(commandArgs(trailingOnly = TRUE)[1])

win <- readRDS(WINDOW_FILE)$windows
expected <- crossing(centre = win$centre[win$used], chr = CHROMOSOMES) %>%
  mutate(stem = sprintf('%sc%02d_chr%d_', SCAN_DIR, centre, chr))
missing <- filter(expected, !file.exists(paste0(stem, 'genes.csv.gz')))
if (nrow(missing)) stop(nrow(missing), ' scan outputs missing, e.g. ', missing$stem[1])
readAll <- \(what) rbindlist(map(paste0(expected$stem, what, '.csv.gz'), \(f) {
  x <- fread(f)
  if (nrow(x)) x[, gene := as.character(gene)]
  x
}))

genes <- as_tibble(readAll('genes'))
nTests <- genes %>% distinct(centre, chr, nTests) %>% count(centre, wt = nTests, name = 'nTestsWindow')
pairs <- as_tibble(readAll('pairs')) %>%
  left_join(distinct(genes, centre, gene, chr), by = c('centre', 'gene'))
t0 <- Sys.time()
called <- callSignals(rename(pairs, set = centre), rename(nTests, set = centre, nTestsSet = nTestsWindow), nCore)   # 00_functions.R
message('BH and plink2 clumping in ', round(as.numeric(difftime(Sys.time(), t0, units = 'mins')), 1), ' min')
sig <- rename(called$sig, centre = set, nTestsWindow = nTestsSet)
fwrite(sig, paste0(CALL_DIR, 'significantPairs.csv.gz'))

genes <- genes %>%
  left_join(count(sig, centre, gene, name = 'nSig'), by = c('centre', 'gene')) %>%
  mutate(nSig = replace_na(nSig, 0L), eGene = nSig > 0) %>%
  left_join(nTests, by = 'centre') %>%
  left_join(select(win, centre, nW = n, nPEER), by = 'centre')
fwrite(genes, paste0(CALL_DIR, 'genes.csv.gz'))
message(nrow(sig), ' significant pairs; eGenes per centre:')
print(count(genes, centre, wt = eGene, name = 'eGenes'), n = Inf)

independent <- rename(called$independent, centre = set) %>%
  left_join(select(genes, centre, gene, nVar, nW, nPEER), by = c('centre', 'gene')) %>%
  select(centre, gene, chr, variant, signal, nSignal, clumpSize, nW, nPEER, nVar,
         beta, se, pNominal, qval)
fwrite(independent, paste0(CALL_DIR, 'independent.csv.gz'))
message(nrow(independent), ' independent signals in ', n_distinct(independent[c('centre', 'gene')]),
        ' centre x eGene combinations')
