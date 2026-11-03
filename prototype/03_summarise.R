# =============================================================================
# 03  Summary: observed hits and vanishing hits against the permuted-DIP null
#
# Totals over windows for each DIP permutation; for each statistic the observed
# value (perm 0), null mean and 95th percentile, and the empirical
# p = (1 + #null >= observed) / (1 + #null).
#
# Run:     Rscript 03_summarise.R            (after every 02 task)
# Outputs: Results_prototype/summary/{overall,permutationTotals}.csv,
#          Results_prototype/summary/nullDistribution.jpeg
# =============================================================================

.a <- commandArgs(FALSE); .f <- sub('^--file=', '', .a[grep('^--file=', .a)])
source(file.path(if (length(.f)) dirname(normalizePath(.f)) else '.', '00_config.R'))

counts <- rbindlist(map(list.files(COUNT_DIR, pattern = '\\.csv$', full.names = TRUE), fread), fill = TRUE)
missing <- setdiff(0:N_PERM, unique(counts$perm))
if (length(missing)) {
  warning(length(missing), ' DIP permutations missing; see missingPermutations.csv')
  fwrite(data.table(perm = missing), paste0(SUMMARY_DIR, 'missingPermutations.csv'))
}
if (0 %in% missing) stop('Observed results (perm 0) are missing')

totals <- counts %>% group_by(perm) %>%
  summarise(across(c(nHits, nVanish), \(x) sum(x, na.rm = TRUE)),
            nWindows = sum(!is.na(nHits)), nBH = sum(qMethod == 'BH', na.rm = TRUE))
overall <- map_dfr(c('nHits', 'nVanish'), \(s) {
  obs <- totals[[s]][totals$perm == 0]; null <- totals[[s]][totals$perm > 0]
  tibble(statistic = s, observed = obs, nullMean = mean(null),
         null95 = unname(quantile(null, 0.95)),
         p = (1 + sum(null >= obs)) / (1 + length(null)), nPerm = length(null))
})
fwrite(totals,  paste0(SUMMARY_DIR, 'permutationTotals.csv'))
fwrite(overall, paste0(SUMMARY_DIR, 'overall.csv'))
print(overall)
if (any(totals$nBH > 0))
  message(sum(totals$nBH), ' windows used BH because the q-value pi0 estimate was 0')

p <- totals %>% filter(perm > 0) %>%
  pivot_longer(c(nHits, nVanish), names_to = 'statistic') %>%
  ggplot(aes(value)) +
  geom_histogram(bins = 30, fill = 'grey70') +
  geom_vline(data = overall, aes(xintercept = observed), colour = 'firebrick') +
  facet_wrap(~statistic, scales = 'free') +
  labs(x = 'Count over all windows', y = 'DIP permutations') +
  theme_classic()
ggsave(paste0(SUMMARY_DIR, 'nullDistribution.jpeg'), p, width = 7, height = 3.5, dpi = 300)
