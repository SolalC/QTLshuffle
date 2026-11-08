# =============================================================================
# 07  Selection-bias control, pass A: reference samples treated as windows
#
# One SLURM task per reference sample k = 1..N_NULL. For every window centre,
# reference k of that centre (nW donors drawn over 24 h, PEER covariates saved
# by 05) goes through exactly the window procedure (same functions as 03/04):
#   1. MatrixEQTL cis scan, chromosome by chromosome (scanChromosome);
#   2. BH over all cis tests of the set, q < PAIR_FDR, and plink2 clumping of
#      each significant gene (callSignals).
# The resulting "pseudo-window" signals are selected like the window signals
# (significance threshold and lead variant) but have no time structure. 08
# tests them against the other references.
#
# Run:     sbatch 07_nullScan.sh       (array 1..N_NULL, after 05)
# Usage:   Rscript 07_nullScan.R <reference sample k> <n cores>
# Outputs: Results/rolling/w<W>h/null/scan/ref<k>.csv.gz
#            one row per pseudo-window signal: centre, ref, gene, chr, variant,
#            signal, nSignal, clumpSize, nVar, beta, se, pNominal, qval
#          Results/rolling/w<W>h/null/scan/ref<k>_summary.csv
#            per centre: nW, cis tests, significant pairs, eGenes, signals
# =============================================================================

.a <- commandArgs(FALSE); .f <- sub('^--file=', '', .a[grep('^--file=', .a)])
source(file.path(if (length(.f)) dirname(normalizePath(.f)) else '.', '00_config.R'))
args <- commandArgs(trailingOnly = TRUE)
k <- as.integer(args[1]); nCore <- as.integer(args[2])
if (k > N_NULL || k > N_REF) stop('Reference sample ', k, ' > N_NULL or N_REF')
outFile <- sprintf('%sref%04d.csv.gz', NULL_SCAN_DIR, k)
sumFile <- sprintf('%sref%04d_summary.csv', NULL_SCAN_DIR, k)
covFile <- sprintf('%sref%04d.rds', REF_COV_DIR, k)
if (!file.exists(covFile)) stop(covFile, ' missing: run 05 for reference ', k, ' first')

tis    <- readRDS(TISSUE_FILE)
refCov <- readRDS(covFile)          # 05: list by centre (centre, donors, cov)
geno   <- openGenotypes(tis$donors)

# 1. Scan every centre's reference k (centres in parallel, chromosomes in turn)
t0 <- Sys.time()
scans <- parallel::mclapply(refCov, \(rc) {
  rows <- match(rc$donors, tis$donors)
  sc <- map(CHROMOSOMES, \(chrom) {
    s <- scanChromosome(tis, geno, rows, rc$cov, chrom)
    list(pairs = mutate(s$pairs, chr = chrom), genes = s$genes, nTests = s$nTests)
  })
  list(centre = rc$centre, nW = length(rows),
       pairs  = bind_rows(map(sc, 'pairs')) %>% mutate(set = rc$centre, .before = 1),
       genes  = bind_rows(map(sc, 'genes')) %>% mutate(set = rc$centre, .before = 1),
       nTests = sum(map_dbl(sc, 'nTests')))
}, mc.cores = nCore, mc.preschedule = FALSE)
failed <- map_lgl(scans, \(x) inherits(x, 'try-error'))
if (any(failed)) stop('Scan failed for ', sum(failed), ' centres: ', as.character(scans[[which(failed)[1]]]))
message(length(scans), ' centres scanned in ', round(as.numeric(difftime(Sys.time(), t0, units = 'mins')), 1), ' min')

# 2. BH and plink2 clumping, as for the windows
pairs  <- bind_rows(map(scans, 'pairs'))
genes  <- bind_rows(map(scans, 'genes'))
nTests <- tibble(set = map_dbl(scans, 'centre'), nTestsSet = map_dbl(scans, 'nTests'))
t0 <- Sys.time()
called <- callSignals(pairs, nTests, nCore)
message('BH and plink2 clumping in ', round(as.numeric(difftime(Sys.time(), t0, units = 'mins')), 1), ' min')

independent <- called$independent %>%
  left_join(select(genes, set, gene, nVar), by = c('set', 'gene')) %>%
  transmute(centre = set, ref = k, gene, chr, variant, signal, nSignal, clumpSize, nVar,
            beta, se, pNominal, qval)
fwrite(independent, outFile)

summary <- tibble(centre = map_dbl(scans, 'centre'), nW = map_int(scans, 'nW'), nTests = nTests$nTestsSet) %>%
  left_join(count(called$sig, centre = set, name = 'sigPairs'), by = 'centre') %>%
  left_join(summarise(group_by(independent, centre), eGenes = n_distinct(gene), signals = n()), by = 'centre') %>%
  mutate(across(c(sigPairs, eGenes, signals), \(v) replace_na(v, 0L)), ref = k, .before = 1)
fwrite(summary, sumFile)
print(as_tibble(summary), n = Inf)
