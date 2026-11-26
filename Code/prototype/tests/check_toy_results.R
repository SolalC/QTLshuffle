# Checks on the prototype toy run.
suppressPackageStartupMessages({ library(tidyverse); library(data.table) })
res <- commandArgs(trailingOnly = TRUE)[1]
fail <- function(...) { message('FAIL: ', ...); quit(status = 1) }

genes <- fread(file.path(res, 'genes.csv'))
if (nrow(genes) != 50) fail('expected 50 genes, found ', nrow(genes))
overall <- fread(file.path(res, 'summary', 'overall.csv'))
if (file.exists(file.path(res, 'summary', 'missingPermutations.csv'))) fail('missing permutations')
v <- overall[statistic == 'nVanish']
if (v$observed <= v$null95) fail('observed vanishing hits not above the null 95th percentile')

# Toy genes 1-10 per chromosome are phase QTL, 11-20 amplitude QTL, 31-40 rhythmic only:
hits <- fread(file.path(res, 'windowHits', 'observed.csv')) %>%
  mutate(k = as.integer(str_sub(gene, 7, 9)),
         arch = case_when(k <= 10 ~ 'phase', k <= 20 ~ 'amplitude', k <= 30 ~ 'mean',
                          k <= 40 ~ 'rhythmic only', TRUE ~ 'null'))
print(hits %>% group_by(arch) %>% summarise(hits = n(), vanishing = sum(vanish)))
message('PASS')
