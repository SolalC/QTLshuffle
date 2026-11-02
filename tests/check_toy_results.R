# Checks on the toy run: completeness, gene exclusion, and that the observed
# vanishing-hit count exceeds the permuted-DIP null where it should.
suppressPackageStartupMessages({ library(tidyverse); library(data.table) })
res <- commandArgs(trailingOnly = TRUE)[1]
fail <- function(...) { message('FAIL: ', ...); quit(status = 1) }

overall <- fread(file.path(res, 'summary', 'overall.csv'))
print(overall)
if (file.exists(file.path(res, 'summary', 'missingPermutations.csv'))) fail('missing permutations')
tisA <- readRDS(file.path(res, 'tissues', 'ToyA.rds'))
if ('ENSG02150.1' %in% tisA$genes$gene_id) fail('CHIRAL gene not dropped')
if (ncol(tisA$cov) != 8) fail('expected 8 fixed covariates, found ', ncol(tisA$cov))

v <- overall[statistic == 'nVanish']
if (v$observed <= v$null95) fail('observed vanishing hits not above the null 95th percentile')

hits <- rbindlist(map(list.files(file.path(res, 'windowHits'), full.names = TRUE), fread)) %>%
  mutate(k = as.integer(str_sub(gene, 7, 9)))
byArch <- hits %>%
  mutate(arch = case_when(k <= 10 ~ 'phase', k <= 20 ~ 'amplitude', k <= 30 ~ 'mean',
                          k <= 40 ~ 'rhythmic only', TRUE ~ 'null')) %>%
  group_by(arch) %>% summarise(hits = n(), vanishing = sum(vanish), .groups = 'drop')
print(byArch)
meanArch <- byArch %>% filter(arch == 'mean')
if (nrow(meanArch) && meanArch$vanishing / meanArch$hits > 0.2) fail('mean eQTL vanish too often')
message('PASS')
