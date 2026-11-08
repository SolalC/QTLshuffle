# =============================================================================
# 07  Associations per window centre
#
# Panel a: on one count scale, per window centre, the number of independent
#          signals (04, BH q < PAIR_FDR, clumped) and the number of those that
#          are time dependent: window beta outside the distribution of the
#          N_REF random 24 h reference samples (06: pEmp < P_EMP, two-sided).
# Panel b: per window, stacked bars summing to 100%: time-dependent signals
#          (pEmp < P_EMP) and time-independent signals (pEmp >= P_EMP).
#          Signals without a reference test (pEmp missing) are left out.
#
# Run:     Rscript 07_plotWindows.R [pEmp threshold, default 0.05] [pEmp | pEmpLarger]
# Outputs: Results/rolling/w<W>h/figures/associationsPerWindow.{pdf,png,csv}
# =============================================================================

.a <- commandArgs(FALSE); .f <- sub('^--file=', '', .a[grep('^--file=', .a)])
source(file.path(if (length(.f)) dirname(normalizePath(.f)) else '.', '00_config.R'))

args  <- commandArgs(TRUE)
P_EMP <- if (length(args) >= 1) as.numeric(args[1]) else 0.05
P_COL <- if (length(args) >= 2) args[2] else 'pEmp'
stopifnot(P_COL %in% c('pEmp', 'pEmpLarger'))

FIG_DIR <- paste0(WIDTH_DIR, 'figures/')
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

x <- fread(paste0(OUT_DIR, 'rollingWindow.csv.gz'), select = c('centre', 'gene', 'variant', 'nRef', P_COL))
setnames(x, P_COL, 'p')
win <- fread(paste0(WIDTH_DIR, 'windows.csv'))

counts <- x %>%
  group_by(centre) %>%
  summarise(signals = n(), eGenes = n_distinct(gene),
            timeDependent = sum(p < P_EMP, na.rm = TRUE), timeIndependent = sum(p >= P_EMP, na.rm = TRUE), untested = sum(is.na(p)),
            minRef = min(nRef), .groups = 'drop') %>%
  right_join(select(filter(win, used), centre, nW = n), by = 'centre') %>%   # skipped windows (n < MIN_WINDOW_N) are left out, not plotted as 0
  mutate(across(c(signals, eGenes, timeDependent, timeIndependent, untested), \(v) replace_na(v, 0L)),
         pctTimeDependent = 100 * timeDependent / (timeDependent + timeIndependent)) %>%
  arrange(centre)
fwrite(counts, paste0(FIG_DIR, 'associationsPerWindow.csv'))
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

stem <- paste0(FIG_DIR, 'associationsPerWindow', if (P_COL != 'pEmp' || P_EMP != 0.05) sprintf('_%s%g', P_COL, P_EMP) else '')
ggsave(paste0(stem, '.pdf'), p, width = 7, height = 7)
ggsave(paste0(stem, '.png'), p, width = 7, height = 7, dpi = 200)
message('Wrote ', stem, '.{pdf,png} and associationsPerWindow.csv')
