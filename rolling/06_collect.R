# =============================================================================
# 06  One table: each independent signal x window centre, with its window
#     association and the N_REF 24 h reference associations
#
# Columns: centre (h), nW, nPEER, gene, variant, signal, window beta/se/p,
# GTEx thresholds and q-value; reference summary: nRef (non-missing),
# refMean, refSD, nSharedMean (window donors in a reference sample, mean),
#   zRef      = (beta - refMean) / refSD
#   pEmp      = (1 + #{|beta_ref - refMean| >= |beta - refMean|}) / (1 + nRef)
#               (two-sided: is the window effect outside the 24 h distribution?)
#   pEmpLarger= (1 + #{beta_ref * sign(beta) >= |beta|}) / (1 + nRef)
#               (one-sided: is the window effect larger than over 24 h?)
# then beta_ref_1 .. beta_ref_<N_REF>. SE and p of every reference
# association are in referenceLong.csv.gz.
#
# Run:     Rscript 06_collect.R        (after every 05 task)
# Outputs: Results_rolling/w<W>h/results/rollingWindow.csv.gz,
#          referenceLong.csv.gz, missingReference.csv (if any)
# =============================================================================

.a <- commandArgs(FALSE); .f <- sub('^--file=', '', .a[grep('^--file=', .a)])
source(file.path(if (length(.f)) dirname(normalizePath(.f)) else '.', '00_config.R'))

signals <- fread(paste0(CALL_DIR, 'independent.csv.gz'))
files <- keep(list.files(REF_DIR, pattern = '\\.csv\\.gz$', full.names = TRUE), \(f) file.size(f) > 0)
ref <- rbindlist(map(files, fread), fill = TRUE)
fwrite(ref, paste0(OUT_DIR, 'referenceLong.csv.gz'))

missing <- distinct(signals, centre) %>% crossing(ref = seq_len(N_REF)) %>%
  anti_join(distinct(ref, centre, ref), by = c('centre', 'ref'))
if (nrow(missing)) {
  warning(nrow(missing), ' centre x reference samples missing; see missingReference.csv')
  fwrite(missing, paste0(OUT_DIR, 'missingReference.csv'))
}

stats <- ref %>%
  left_join(select(signals, centre, gene, variant, betaW = beta), by = c('centre', 'gene', 'variant')) %>%
  group_by(centre, gene, variant) %>%
  summarise(nRef = sum(!is.na(beta)), refMean = mean(beta, na.rm = TRUE), refSD = sd(beta, na.rm = TRUE),
            nSharedMean = mean(nShared),
            pEmp = (1 + sum(abs(beta - refMean) >= abs(betaW[1] - refMean), na.rm = TRUE)) / (1 + nRef),
            pEmpLarger = (1 + sum(beta * sign(betaW[1]) >= abs(betaW[1]), na.rm = TRUE)) / (1 + nRef),
            .groups = 'drop')
wide <- ref %>% select(centre, gene, variant, ref, beta) %>%
  mutate(ref = paste0('beta_ref_', ref)) %>%
  pivot_wider(names_from = ref, values_from = beta)
wide <- wide[, c('centre', 'gene', 'variant', intersect(paste0('beta_ref_', seq_len(N_REF)), names(wide)))]

out <- signals %>%
  left_join(stats, by = c('centre', 'gene', 'variant')) %>%
  mutate(zRef = (beta - refMean) / refSD, .after = nSharedMean) %>%
  left_join(wide, by = c('centre', 'gene', 'variant')) %>%
  arrange(gene, centre, signal)
fwrite(out, paste0(OUT_DIR, 'rollingWindow.csv.gz'))
message(nrow(out), ' signal x centre rows; ', sum(out$pEmp < 0.05, na.rm = TRUE), ' with pEmp < 0.05')
