# =============================================================================
# 08  Example signals: window beta against its N_REF random 24 h references
#
# For one time-dependent and one time-independent signal, the histogram of
# the N_REF reference betas (same gene, variant and n, donors drawn over 24 h)
# with a vertical line at the window beta.
#   time dependent  : smallest pEmp (ties: largest |zRef|)
#   time independent: largest pEmp  (ties: smallest |zRef|)
# Both are significant in their window (04), so the time-independent example
# is a window association whose beta is typical of the 24 h distribution.
# A signal can be chosen instead: Rscript 08_plotExamples.R <centre> <gene> <variant>
# (the row is then plotted on its own, whatever its pEmp).
#
# Run:     Rscript 08_plotExamples.R [centre gene variant]
# Outputs: Results/rolling/w<W>h/figures/exampleSignals.{pdf,png,csv}
#          (or exampleSignal_c<centre>_<gene>_<variant>.* for a chosen signal)
# =============================================================================

.a <- commandArgs(FALSE); .f <- sub('^--file=', '', .a[grep('^--file=', .a)])
source(file.path(if (length(.f)) dirname(normalizePath(.f)) else '.', '00_config.R'))

P_EMP   <- 0.05
FIG_DIR <- paste0(WIDTH_DIR, 'figures/')
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
args <- commandArgs(TRUE)

x <- fread(paste0(OUT_DIR, 'rollingWindow.csv.gz'))
x <- x[!is.na(pEmp)]

if (length(args) == 3) {
  ex <- x[centre == as.numeric(args[1]) & gene == args[2] & variant == args[3]]
  if (!nrow(ex)) stop('No tested signal ', paste(args, collapse = ' '))
  ex[, example := if_else(pEmp < P_EMP, 'Time dependent', 'Time independent')]
  stem <- sprintf('%sexampleSignal_c%02d_%s_%s', FIG_DIR, ex$centre, ex$gene, gsub(':', '_', ex$variant))
} else {
  ex <- bind_rows(
    x[order(pEmp, -abs(zRef))][1][, example := 'Time dependent'],
    x[order(-pEmp, abs(zRef))][1][, example := 'Time independent'])
  stem <- paste0(FIG_DIR, 'exampleSignals')
}

info <- ex %>%
  select(example, centre, gene, variant, nW, beta, se, pNominal, qval, nRef, refMean, refSD, zRef, pEmp, pEmpLarger) %>%
  mutate(label = sprintf('%s\n%s, %s, window %02d h (n = %d)\nwindow beta = %.3f; reference mean %.3f, SD %.3f; z = %.2f; pEmp = %.3g',
                         example, gene, variant, centre, nW, beta, refMean, refSD, zRef, pEmp))
fwrite(select(info, -label), paste0(stem, '.csv'))
print(as_tibble(select(info, -label)), width = Inf)

refs <- ex %>%
  select(example, starts_with('beta_ref_')) %>%
  pivot_longer(-example, values_to = 'betaRef', values_drop_na = TRUE) %>%
  left_join(select(info, example, label), by = 'example')

p <- ggplot(refs, aes(betaRef)) +
  geom_histogram(bins = 50, fill = '#b9c3cc', colour = 'white', linewidth = 0.2) +
  geom_vline(data = info, aes(xintercept = beta), colour = '#b2443a', linewidth = 0.9) +
  facet_wrap(~ label, ncol = 1, scales = 'free') +
  labs(x = 'Effect size (beta)', y = sprintf('Reference samples (of %d)', N_REF),
       title = sprintf('%s: %g h window vs %d random 24 h samples', TISSUE, WINDOW_WIDTH_H, N_REF),
       subtitle = 'Grey: reference betas (same n, donors drawn over 24 h). Red: window beta.') +
  theme_bw(base_size = 11) +
  theme(strip.background = element_blank(), strip.text = element_text(hjust = 0),
        panel.grid.minor = element_blank())

ggsave(paste0(stem, '.pdf'), p, width = 7, height = 3.2 * nrow(info) + 0.8)
ggsave(paste0(stem, '.png'), p, width = 7, height = 3.2 * nrow(info) + 0.8, dpi = 200)
message('Wrote ', stem, '.{pdf,png,csv}')
