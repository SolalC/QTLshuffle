# =============================================================================
# 04  Summary: observed vanishing hits against the permuted-DIP null
#
# Totals over the selected tissues and windows for each permutation, then for
# each statistic: observed (perm 0), null mean and 95th percentile, and the
# empirical p = (1 + #null >= observed) / (1 + #null). Also per tissue.
# Only permutations complete for every tissue enter the totals.
#
# Run:     Rscript 04_summarise.R            (after every 03 task)
# Outputs: Results/main/summary/{overall,perTissue,permutationTotals}.csv,
#          Results/main/summary/nullDistribution.jpeg
# =============================================================================

.a <- commandArgs(FALSE); .f <- sub('^--file=', '', .a[grep('^--file=', .a)])
source(file.path(if (length(.f)) dirname(normalizePath(.f)) else '.', '00_config.R'))

tissues <- fread(paste0(TISSUE_DIR, 'selectedTissues.csv'))
counts  <- rbindlist(map(list.files(COUNT_DIR, pattern = '\\.csv$', full.names = TRUE), fread))

expected <- crossing(tissue = tissues$tissue, perm = 0:N_PERM)
present  <- counts %>% distinct(tissue, perm)
missing  <- anti_join(expected, present, by = c('tissue', 'perm'))
if (nrow(missing)) {
  warning(nrow(missing), ' tissue x permutation results missing; see missingPermutations.csv')
  fwrite(missing, paste0(SUMMARY_DIR, 'missingPermutations.csv'))
}
if (0 %in% missing$perm) stop('Observed results (perm 0) are missing for some tissues')
completePerms <- setdiff(0:N_PERM, missing$perm)

stats <- c('nHits', 'nVanish', 'nVanishRhythmic', 'nVanishNonRhythmic')
addNonRhythmic <- function(d) d %>% mutate(nVanishNonRhythmic = nVanish - nVanishRhythmic)

empirical <- function(d) {
  map_dfr(stats, \(s) {
    obs  <- d[[s]][d$perm == 0]
    null <- d[[s]][d$perm > 0]
    tibble(statistic = s, observed = obs, nullMean = mean(null),
           null95 = unname(quantile(null, 0.95)),
           p = (1 + sum(null >= obs)) / (1 + length(null)), nPerm = length(null))
  })
}

totals <- counts %>% filter(perm %in% completePerms) %>%
  group_by(perm) %>%
  summarise(across(c(nHits, nVanish, nVanishRhythmic), \(x) sum(x, na.rm = TRUE))) %>%
  addNonRhythmic()
perTissue <- counts %>% filter(perm %in% completePerms) %>%
  group_by(tissue, perm) %>%
  summarise(across(c(nHits, nVanish, nVanishRhythmic), \(x) sum(x, na.rm = TRUE)),
            .groups = 'drop') %>%
  addNonRhythmic() %>%
  group_by(tissue) %>% group_modify(\(d, k) empirical(d)) %>% ungroup()

overall <- empirical(totals)
fwrite(totals,    paste0(SUMMARY_DIR, 'permutationTotals.csv'))
fwrite(overall,   paste0(SUMMARY_DIR, 'overall.csv'))
fwrite(perTissue, paste0(SUMMARY_DIR, 'perTissue.csv'))
print(overall)

p <- ggplot(filter(totals, perm > 0), aes(nVanish)) +
  geom_histogram(bins = 40, fill = 'grey70') +
  geom_vline(xintercept = totals$nVanish[totals$perm == 0], colour = 'firebrick') +
  labs(x = 'Vanishing hits (all tissues and windows)', y = 'Permutations',
       subtitle = sprintf('Observed = %d; empirical p = %.3g',
                          totals$nVanish[totals$perm == 0],
                          overall$p[overall$statistic == 'nVanish'])) +
  theme_classic()
ggsave(paste0(SUMMARY_DIR, 'nullDistribution.jpeg'), p, width = 5, height = 3.5, dpi = 300)
