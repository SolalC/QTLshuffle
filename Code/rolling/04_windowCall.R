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
library(parallel)
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
  left_join(nTests, by = 'centre') %>%
  group_by(centre) %>% mutate(qval = p.adjust(pNominal, 'BH', n = nTestsWindow[1])) %>% ungroup()
sig <- filter(pairs, qval < PAIR_FDR) %>%
  left_join(distinct(genes, centre, gene, chr), by = c('centre', 'gene')) %>%
  arrange(centre, gene, pNominal)
maxP <- sig %>% group_by(centre) %>% summarise(maxP = max(pNominal))
if (any(maxP$maxP > PAIR_P_STORE / 10))
  warning('Significant pairs reach p = ', signif(max(maxP$maxP), 2), ', close to PAIR_P_STORE (',
          PAIR_P_STORE, '); raise PAIR_P_STORE in 00_config.R and rerun 03')
fwrite(sig, paste0(CALL_DIR, 'significantPairs.csv.gz'))

genes <- genes %>%
  left_join(count(sig, centre, gene, name = 'nSig'), by = c('centre', 'gene')) %>%
  mutate(nSig = replace_na(nSig, 0L), eGene = nSig > 0) %>%
  left_join(nTests, by = 'centre') %>%
  left_join(select(win, centre, nW = n, nPEER), by = 'centre')
fwrite(genes, paste0(CALL_DIR, 'genes.csv.gz'))
message(nrow(sig), ' significant pairs; eGenes per centre:')
print(count(genes, centre, wt = eGene, name = 'eGenes'), n = Inf)

# plink2 --clump for one centre x eGene; returns index variants and clump sizes.
clump <- function(x) {
  dir <- tempfile('clump_'); dir.create(dir); on.exit(unlink(dir, recursive = TRUE))
  fwrite(tibble(ID = x$variant, P = pmax(x$pNominal, 1e-300)), file.path(dir, 'assoc.tsv'), sep = '\t')
  status <- system2(PLINK2, c('--bfile', PLINK_CHR(x$chr[1]), '--clump', file.path(dir, 'assoc.tsv'),
                              '--clump-p1', 1, '--clump-p2', 1, '--clump-r2', CLUMP_R2,
                              '--clump-kb', CLUMP_KB, '--threads', 1, '--memory', 4000,
                              '--out', file.path(dir, 'out')), stdout = FALSE, stderr = FALSE)
  out <- file.path(dir, 'out.clumps')
  if (status != 0 || !file.exists(out)) stop('plink2 --clump failed for ', x$centre[1], ' h, ', x$gene[1])
  fread(out) %>% transmute(variant = ID, clumpSize = TOTAL + 1L)
}
groups <- split(sig, interaction(sig$centre, sig$gene, drop = TRUE))
t0 <- Sys.time()
index <- mclapply(groups, \(x) mutate(clump(x), centre = x$centre[1], gene = x$gene[1]),
                  mc.cores = nCore, mc.preschedule = FALSE)
failed <- map_lgl(index, \(x) inherits(x, 'try-error'))
if (any(failed)) stop(sum(failed), ' clumping runs failed, e.g. ', as.character(index[[which(failed)[1]]]))
message(length(groups), ' plink2 clumping runs in ', round(as.numeric(difftime(Sys.time(), t0, units = 'mins')), 1), ' min')

independent <- bind_rows(index) %>%
  inner_join(sig, by = c('centre', 'gene', 'variant')) %>%
  arrange(centre, gene, pNominal) %>%
  group_by(centre, gene) %>% mutate(signal = row_number(), nSignal = n()) %>% ungroup() %>%
  left_join(select(genes, centre, gene, nVar, nW, nPEER), by = c('centre', 'gene')) %>%
  select(centre, gene, chr, variant, signal, nSignal, clumpSize, nW, nPEER, nVar,
         beta, se, pNominal, qval)
fwrite(independent, paste0(CALL_DIR, 'independent.csv.gz'))
message(nrow(independent), ' independent signals in ', n_distinct(independent[c('centre', 'gene')]),
        ' centre x eGene combinations')
