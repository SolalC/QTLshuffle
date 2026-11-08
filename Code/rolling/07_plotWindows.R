# =============================================================================
# 07  Figure: example signals (left) and associations per window centre (right)
#
#   a  time-independent example: histogram of the N_REF reference betas
#      (same gene, variant and n, donors drawn over 24 h), red line at the
#      window beta. Default: largest pEmp (ties: smallest |zRef|).
#   b  per window centre, on one count scale: independent signals (04, BH
#      q < PAIR_FDR, clumped) and those that are time dependent (window beta
#      outside the 24 h reference distribution, 06: pEmp < P_EMP, two-sided).
#   c  time-dependent example, as a. Default: smallest pEmp (ties: largest |zRef|).
#   d  per window, stacked bars summing to 100%: time dependent (pEmp < P_EMP)
#      and time independent (pEmp >= P_EMP). Signals without a reference test
#      (pEmp missing) are left out.
# Both examples are significant in their window (04).
#
# Run:     Rscript 07_plotWindows.R [p=0.05] [test=pEmp|pEmpLarger]
#            [independent=centre,gene,variant] [dependent=centre,gene,variant]
#          (independent= / dependent= replace the automatic examples)
# Outputs: Results/rolling/w<W>h/figures/rollingWindow.{pdf,png},
#          associationsPerWindow.csv, exampleSignals.csv
#          (suffix _<test><p> when the test or threshold is not the default)
# =============================================================================

.a <- commandArgs(FALSE); .f <- sub('^--file=', '', .a[grep('^--file=', .a)])
source(file.path(if (length(.f)) dirname(normalizePath(.f)) else '.', '00_config.R'))
library(patchwork)

args   <- commandArgs(TRUE)
argVal <- function(key, default) {
  v <- sub(paste0('^', key, '='), '', args[startsWith(args, paste0(key, '='))])
  if (length(v)) v[1] else default
}
P_EMP  <- as.numeric(argVal('p', 0.05))
P_COL  <- argVal('test', 'pEmp')
EX_IND <- argVal('independent', NA)
EX_DEP <- argVal('dependent', NA)
stopifnot(P_COL %in% c('pEmp', 'pEmpLarger'), P_EMP > 0, P_EMP < 1)
SUFFIX <- if (P_COL != 'pEmp' || P_EMP != 0.05) sprintf('_%s%g', P_COL, P_EMP) else ''

FIG_DIR <- paste0(WIDTH_DIR, 'figures/')
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

full <- fread(paste0(OUT_DIR, 'rollingWindow.csv.gz'))
full[, p := get(P_COL)]
win  <- fread(paste0(WIDTH_DIR, 'windows.csv'))

LAB_ALL <- sprintf('All signals (window q < %g)', PAIR_FDR)
LAB_DEP <- 'Time dependent'
LAB_IND <- 'Time independent'
SERIES  <- c(LAB_ALL, LAB_DEP, LAB_IND)
COLOURS <- setNames(c('#2b5d8a', '#b2443a', '#b9c3cc'), SERIES)
X_TIME  <- scale_x_continuous(breaks = seq(0, 21, 3), limits = c(-0.5, 23.5), expand = expansion(0))
THEME   <- theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank())

# --- Counts per window (b, d) ----------------------------------------------------------
counts <- full %>%
  group_by(centre) %>%
  summarise(signals = n(), eGenes = n_distinct(gene),
            timeDependent = sum(p < P_EMP, na.rm = TRUE), timeIndependent = sum(p >= P_EMP, na.rm = TRUE),
            untested = sum(is.na(p)), minRef = min(nRef), .groups = 'drop') %>%
  right_join(select(filter(win, used), centre, nW = n), by = 'centre') %>%   # skipped windows (n < MIN_WINDOW_N) are left out, not plotted as 0
  mutate(across(c(signals, eGenes, timeDependent, timeIndependent, untested), \(v) replace_na(v, 0L)),
         pctTimeDependent = 100 * timeDependent / (timeDependent + timeIndependent)) %>%
  arrange(centre)
fwrite(counts, paste0(FIG_DIR, 'associationsPerWindow', SUFFIX, '.csv'))
print(as_tibble(counts), n = Inf)

lines <- bind_rows(
  transmute(counts, centre, nW, series = LAB_ALL, value = signals),
  transmute(counts, centre, nW, series = LAB_DEP, value = timeDependent)) %>%
  mutate(series = factor(series, SERIES))
bars <- bind_rows(
  transmute(counts, centre, series = LAB_DEP, value = timeDependent),
  transmute(counts, centre, series = LAB_IND, value = timeIndependent)) %>%
  group_by(centre) %>% mutate(value = 100 * value / sum(value)) %>% ungroup() %>%
  filter(!is.na(value)) %>%
  mutate(series = factor(series, rev(SERIES)))   # time dependent at the bottom

pCount <- ggplot(lines, aes(centre, value, colour = series)) +
  geom_line() + geom_point(aes(size = nW)) +
  X_TIME + scale_y_continuous(limits = c(0, NA)) +
  scale_colour_manual(values = COLOURS, breaks = c(LAB_ALL, LAB_DEP), name = NULL) +
  scale_size_continuous(range = c(0.8, 3), name = 'Window donors') +
  labs(x = 'Window centre (h, DIP time)', y = 'Number of associations') + THEME
pBar <- ggplot(bars, aes(centre, value, fill = series)) +
  geom_col(width = 0.85) +
  X_TIME + scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25), expand = expansion(0)) +
  scale_fill_manual(values = COLOURS, breaks = c(LAB_DEP, LAB_IND), name = NULL) +
  labs(x = 'Window centre (h, DIP time)', y = 'Associations (%)') + THEME

# --- Examples (a, c) -------------------------------------------------------------------
tested <- full[!is.na(p)]
pick <- function(key, auto, label) {
  ex <- if (is.na(key)) auto else {
    k <- strsplit(key, ',')[[1]]
    tested[centre == as.numeric(k[1]) & gene == k[2] & variant == k[3]]
  }
  if (nrow(ex) != 1) stop('No tested signal ', key)
  ex[, example := label]
}
ex <- bind_rows(
  pick(EX_IND, tested[order(-p, abs(zRef))][1], LAB_IND),
  pick(EX_DEP, tested[order(p, -abs(zRef))][1], LAB_DEP))
info <- ex %>%
  select(example, centre, gene, variant, nW, beta, se, pNominal, qval, nRef, refMean, refSD, zRef, pEmp, pEmpLarger, p) %>%
  mutate(label = sprintf('%s: %s, %s\nwindow %02d h (n = %d); z = %.2f; %s = %.3g',
                         example, gene, variant, centre, nW, round(zRef, 2) + 0, P_COL, p))
fwrite(select(info, -label, -p), paste0(FIG_DIR, 'exampleSignals', SUFFIX, '.csv'))
print(as_tibble(select(info, -label)), width = Inf)
if (with(info, any((example == LAB_DEP) != (p < P_EMP))))
  warning('A chosen example is not classified as its panel says at ', P_COL, ' < ', P_EMP)

histogram <- function(lab) {
  i    <- filter(info, example == lab)
  refs <- ex[example == lab] %>% select(starts_with('beta_ref_')) %>%
    pivot_longer(everything(), values_to = 'betaRef', values_drop_na = TRUE)
  ggplot(refs, aes(betaRef)) +
    geom_histogram(bins = 50, fill = COLOURS[LAB_IND], colour = 'white', linewidth = 0.2) +
    geom_vline(xintercept = i$beta, colour = COLOURS[LAB_DEP], linewidth = 0.9) +
    labs(x = 'Effect size (beta); red: window', y = sprintf('Reference samples (of %d)', N_REF),
         subtitle = i$label) +
    THEME + theme(plot.subtitle = element_text(size = 8))
}

# --- Figure ---------------------------------------------------------------------------
fig <- (histogram(LAB_IND) | pCount) / (histogram(LAB_DEP) | pBar) +
  plot_layout(widths = c(1, 1.4)) +
  plot_annotation(tag_levels = 'a') &
  theme(legend.position = 'right', legend.justification = 'left', plot.tag = element_text(face = 'bold'))

stem <- paste0(FIG_DIR, 'rollingWindow', SUFFIX)
ggsave(paste0(stem, '.pdf'), fig, width = 12, height = 7)
ggsave(paste0(stem, '.png'), fig, width = 12, height = 7, dpi = 200)
message('Wrote ', stem, '.{pdf,png}, associationsPerWindow', SUFFIX, '.csv, exampleSignals', SUFFIX, '.csv')
