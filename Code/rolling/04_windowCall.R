# =============================================================================
# 04  eGenes, significant pairs and independent signals per window centre
#
# For each centre, across all genes of all chromosomes (GTEx protocol, tensorQTL
# post.calculate_qvalues):
#   1. Storey q-value of pBeta (lambda = QVALUE_LAMBDA); eGene: q <= GENE_FDR;
#   2. gene-level nominal threshold qbeta(pthreshold, shape1, shape2);
#      significant pairs: cis variants of an eGene with nominal p <= threshold;
#   3. independent signals: greedy clumping of each eGene's significant variants
#      at r2 > CLUMP_R2 (as plink --clump, no distance limit within the cis
#      window). LD (r2) is computed on all the tissue's donors, not the window
#      donors: with ~90 donors, r2 between unlinked variants is ~1/90 on
#      average, above the 0.01 cut-off, so window-donor LD would merge
#      independent signals at random.
#
# Run:     sbatch 04_windowCall.sh      (after every 03 task)
# Usage:   Rscript 04_windowCall.R <n cores>
# Outputs: Results/rolling/w<W>h/calls/
#            genes.csv.gz            every gene tested, per centre, with q-value
#            significantPairs.csv.gz significant pairs with their index variant
#            independent.csv.gz      index variants (one row per signal)
#          GENO_DIR/<tissue>/w<W>h_indexGenotypes.rds: raw genotypes of the index
#            variants (all donors), for 05
# =============================================================================

.a <- commandArgs(FALSE); .f <- sub('^--file=', '', .a[grep('^--file=', .a)])
source(file.path(if (length(.f)) dirname(normalizePath(.f)) else '.', '00_config.R'))
library(parallel)
nCore <- as.integer(commandArgs(trailingOnly = TRUE)[1])

win <- readRDS(WINDOW_FILE)$windows
used <- filter(win, used)
expected <- crossing(centre = used$centre, chr = CHROMOSOMES) %>%
  mutate(file = sprintf('%sc%02d_chr%d_genes.csv.gz', SCAN_DIR, centre, chr))
missing <- filter(expected, !file.exists(file))
if (nrow(missing)) stop(nrow(missing), ' scan outputs missing, e.g. ', missing$file[1])

genes <- rbindlist(map(expected$file, fread)) %>%
  group_by(centre) %>% group_modify(\(g, k) gtexQvalue(g)) %>% ungroup() %>%
  left_join(select(win, centre, nW = n, nPEER), by = 'centre')
fwrite(genes, paste0(CALL_DIR, 'genes.csv.gz'))
eGenes <- filter(genes, eGene)
message('eGenes per centre:'); print(count(genes, centre, wt = eGene, name = 'eGenes'), n = Inf)
over <- filter(eGenes, pvalNominalThreshold >= PAIR_P_STORE)
if (nrow(over)) warning(nrow(over), ' eGenes have a nominal threshold >= PAIR_P_STORE (',
                        PAIR_P_STORE, '); their significant pairs are truncated')

pairs <- rbindlist(map(sprintf('%sc%02d_chr%d_pairs.csv.gz', SCAN_DIR, expected$centre, expected$chr), \(f) {
  x <- fread(f)
  if (nrow(x)) x[, gene := as.character(gene)]
  x
})) %>%
  inner_join(select(eGenes, centre, gene, chr, pvalNominalThreshold), by = c('centre', 'gene')) %>%
  filter(pNominal <= pvalNominalThreshold) %>%
  arrange(centre, gene, pNominal)
message(nrow(pairs), ' significant pairs')

# Clumping, one chromosome at a time (genotypes of all the tissue's donors).
clumped <- map_dfr(CHROMOSOMES, \(chrom) {
  pc <- filter(pairs, chr == chrom)
  if (!nrow(pc)) return(NULL)
  G <- readRDS(GENO_FILE(chrom))$G
  groups <- split(pc, interaction(pc$centre, pc$gene, drop = TRUE))
  res <- mclapply(groups, \(x) mutate(x, indexVariant = clumpGene(dosage(G[, x$variant, drop = FALSE]))),
                  mc.cores = nCore)
  failed <- map_lgl(res, \(x) inherits(x, 'try-error'))
  if (any(failed)) stop('Clumping failed: ', as.character(res[[which(failed)[1]]]))
  bind_rows(res)
})
if (!nrow(clumped)) clumped <- mutate(pairs, indexVariant = character())
fwrite(clumped, paste0(CALL_DIR, 'significantPairs.csv.gz'))

independent <- clumped %>%
  group_by(centre, gene, indexVariant) %>% mutate(clumpSize = n()) %>% ungroup() %>%
  filter(variant == indexVariant) %>%
  group_by(centre, gene) %>% mutate(signal = row_number(), nSignal = n()) %>% ungroup() %>%
  left_join(select(genes, centre, gene, qval, pBeta, nVar, nW, nPEER), by = c('centre', 'gene')) %>%
  select(centre, gene, chr, variant, signal, nSignal, clumpSize, nW, nPEER, nVar,
         beta, se, pNominal, pvalNominalThreshold, pBeta, qval)
fwrite(independent, paste0(CALL_DIR, 'independent.csv.gz'))
message(nrow(independent), ' independent signals in ', n_distinct(independent[c('centre', 'gene')]),
        ' centre x eGene combinations')

indexG <- map(unique(independent$chr), \(chrom) {
  G <- readRDS(GENO_FILE(chrom))$G
  G[, unique(independent$variant[independent$chr == chrom]), drop = FALSE]
})
saveRDS(if (length(indexG)) do.call(cbind, indexG) else NULL, INDEX_GENO_FILE)
