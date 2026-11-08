# =============================================================================
# 07  Figures: associations per window centre, and example signals
#
# Figure 1 (associationsPerWindow):
#   panel a: on one count scale, per window centre, the number of independent
#            signals (04, BH q < PAIR_FDR, clumped) and the number of those
#            that are time dependent: window beta outside the distribution of
#            the N_REF random 24 h reference samples (06: pEmp < P_EMP,
#            two-sided).
#   panel b: per window, stacked bars summing to 100%: time-dependent signals
#            (pEmp < P_EMP) and time-independent signals (pEmp >= P_EMP).
#            Signals without a reference test (pEmp missing) are left out.
# Figure 2 (exampleSignals): for one time-dependent and one time-independent
#   signal, the histogram of the N_REF reference betas (same gene, variant
#   and n, donors drawn over 24 h) with a vertical line at the window beta.
#     time dependent  : smallest pEmp (ties: largest |zRef|)
#     time independent: largest pEmp  (ties: smallest |zRef|)
#   Both are significant in their window (04). With example=<centre>,<gene>,<variant>
#   that signal is plotted instead (exampleSignal_c<centre>_<gene>_<variant>.*).
#
# Run:     Rscript 07_plotWindows.R [p=0.05] [test=pEmp|pEmpLarger] [example=centre,gene,variant]
# Outputs: Results/rolling/w<W>h/figures/associationsPerWindow.{pdf,png,csv},
#          exampleSignals.{pdf,png,csv}
# =============================================================================

.a <- commandArgs(FALSE); .f <- sub('^--file=', '', .a[grep('^--file=', .a)])
source(file.path(if (length(.f)) dirname(normalizePath(.f)) else '.', '00_config.R'))

args    <- commandArgs(TRUE)
argVal  <- function(key, default) {
  v <- sub(paste0('^', key, '='), '', args[startsWith(args, paste0(key, '='))])
  if (length(v)) v[1] else default
}
P_EMP   <- as.numeric(argVal('p', 0.05))
P_COL   <- argVal('test', 'pEmp')
EXAMPLE <- argVal('example', NA)
stopifnot(P_COL %in% c('pEmp', 'pEmpLarger'), P_EMP > 0, P_EMP < 1)
SUFFIX  <- if (P_COL != 'pEmp' || P_EMP != 0.05) sprintf('_%s%g', P_COL, P_EMP) else ''   # non-default test: separate files

FIG_DIR <- paste0(WIDTH_DIR, 'figures/')
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

full <- fread(paste0(OUT_DIR, 'rollingWindow.csv.gz'))
x    <- full[, .(centre, gene, variant, nRef, p = get(P_COL))]
win  <- fread(paste0(WIDTH_DIR, 'windows.csv'))

# --- Figure 1: associations per window ------------------------------------------
counts <- x %>%
  group_by(centre) %>%
  summarise(signals = n(), eGenes = n_distinct(gene),
            timeDependent = sum(p < P_EMP, na.rm = TRUE), timeIndependent = sum(p >= P_EMP, na.rm = TRUE), untested = sum(is.na(p)),
            minRef = min(nRef), .groups = 'drop') %>%
  right_join(select(filter(win, used), centre, nW = n), by = 'centre') %>%   # skipped windows (n < MIN_WINDOW_N) are left out, not plotted as 0
  mutate(across(c(signals, eGenes, timeDependent, timeIndependent, untested), \(v) replace_na(v, 0L)),
         pctTimeDependent = 100 * timeDependent / (timeDependent + timeIndependent)) %>%
  arrange(centre)
fwrite(counts, paste0(FIG_DIR, 'associationsPerWindow', SUFFIX, '.csv'))
print(as_tibble(counts), n = Inf)

LAB_ALL <- sprintf('All signals (window q < %g)', PAIR_FDR)
LAB_DEP <- 'Time dependent'
LAB_IND <- 'Time independent'
PANEL_A <- 'a  Number of associations'
PANEL_B <- sprintf('b  Proportion of associations (%s threshold %g)', P_COL, P_EMP)
PANELS  <- c(PANEL_A, PANEL_B)
SERIES  <- c(LAB_ALL, LAB_DEP, LAB_IND)

lines <- bind_rows(
  transmute(counts, centre, nW, series = LAB_ALL, value = signals),
  transmute(counts, centre, nW, series = LAB_DEP, value = timeDependent)) %>%
  mutate(panel = factor(PANEL_A, PANELS), series = factor(series, SERIES))
bars <- bind_rows(
  transmute(counts, centre, series = LAB_DEP, value = timeDependent),
  transmute(counts, centre, series = LAB_IND, value = timeIndependent)) %>%
  group_by(centre) %>% mutate(value = 100 * value / sum(value)) %>% ungroup() %>%
  filter(!is.na(value)) %>%
  mutate(panel = factor(PANEL_B, PANELS), series = factor(series, rev(SERIES)))   # time dependent at the bottom

p <- ggplot(mapping = aes(centre, value)) +
  geom_col(data = bars, aes(fill = series), width = 0.85) +
  geom_line(data = lines, aes(colour = series)) +
  geom_point(data = lines, aes(colour = series, size = nW)) +
  facet_wrap(~ panel, ncol = 1, scales = 'free_y') +
  scale_x_continuous(breaks = seq(0, 21, 3), limits = c(-0.5, 23.5), expand = expansion(0)) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  scale_colour_manual(values = setNames(c('#2b5d8a', '#b2443a', '#b9c3cc'), SERIES), breaks = SERIES,
                      aesthetics = c('colour', 'fill'), name = NULL) +
  scale_size_continuous(range = c(1, 3.5), name = 'Window donors') +
  labs(x = 'Window centre (h, DIP time)', y = 'Number (a) / % (b)',
       title = sprintf('%s: %g h rolling windows', TISSUE, WINDOW_WIDTH_H)) +
  theme_bw(base_size = 11) +
  theme(strip.background = element_blank(), strip.text = element_text(hjust = 0, face = 'bold'),
        panel.grid.minor = element_blank(), legend.position = 'bottom', legend.box = 'vertical')

stem <- paste0(FIG_DIR, 'associationsPerWindow', SUFFIX)
ggsave(paste0(stem, '.pdf'), p, width = 7, height = 7)
ggsave(paste0(stem, '.png'), p, width = 7, height = 7, dpi = 200)
message('Wrote ', stem, '.{pdf,png,csv}')

# --- Figure 2: example signals --------------------------------------------------
tested <- full[!is.na(get(P_COL))][, p := get(P_COL)]
if (!is.na(EXAMPLE)) {
  key <- strsplit(EXAMPLE, ',')[[1]]
  ex  <- tested[centre == as.numeric(key[1]) & gene == key[2] & variant == key[3]]
  if (!nrow(ex)) stop('No tested signal ', EXAMPLE)
  ex[, example := if_else(p < P_EMP, 'Time dependent', 'Time independent')]
  stem <- sprintf('%sexampleSignal_c%02d_%s_%s%s', FIG_DIR, ex$centre, ex$gene, gsub(':', '_', ex$variant), SUFFIX)
} else {
  ex <- bind_rows(
    tested[order(p, -abs(zRef))][1][, example := 'Time dependent'],
    tested[order(-p, abs(zRef))][1][, example := 'Time independent'])
  stem <- paste0(FIG_DIR, 'exampleSignals', SUFFIX)
}

info <- ex %>%
  select(example, centre, gene, variant, nW, beta, se, pNominal, qval, nRef, refMean, refSD, zRef, pEmp, pEmpLarger, p) %>%
  mutate(label = sprintf('%s\n%s, %s, window %02d h (n = %d)\nwindow beta = %.3f; reference mean %.3f, SD %.3f; z = %.2f; %s = %.3g',
                         example, gene, variant, centre, nW, beta, refMean, refSD, zRef, P_COL, p))
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
