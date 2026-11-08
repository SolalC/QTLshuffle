# =============================================================================
# 07  Associations per window centre
#
# Panel A: on one count scale, per window centre, the number of independent
#          signals (04, BH q < PAIR_FDR, clumped) and the number of those
#          whose window beta lies outside the distribution of the N_REF random
#          24 h reference samples (06: pEmp < P_EMP, two-sided).
# Panel B: the second count as a percentage of the first. The dashed line is
#          the percentage expected if no signal depended on time (100 x P_EMP).
#          Selection on the window beta (winner's curse) inflates panel B
#          above that line even without any time effect, so the line is a
#          floor, not a null.
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
            timeSpecific = sum(p < P_EMP, na.rm = TRUE), untested = sum(is.na(p)),
            minRef = min(nRef), .groups = 'drop') %>%
  right_join(select(filter(win, used), centre, nW = n), by = 'centre') %>%   # skipped windows (n < MIN_WINDOW_N) are left out, not plotted as 0
  mutate(across(c(signals, eGenes, timeSpecific, untested), \(v) replace_na(v, 0L)),
         pctTime = 100 * timeSpecific / (signals - untested)) %>%
  arrange(centre)
fwrite(counts, paste0(FIG_DIR, 'associationsPerWindow.csv'))
print(as_tibble(counts), n = Inf)

LAB_ALL  <- sprintf('Independent signals (window q < %g)', PAIR_FDR)
LAB_TIME <- sprintf('Outside 24 h reference (%s < %g)', P_COL, P_EMP)
PANEL_A  <- 'A  Number of associations'
PANEL_B  <- sprintf('B  %% of signals outside 24 h reference (%s < %g)', P_COL, P_EMP)
long <- bind_rows(
  transmute(counts, centre, nW, panel = PANEL_A, series = LAB_ALL,  value = signals),
  transmute(counts, centre, nW, panel = PANEL_A, series = LAB_TIME, value = timeSpecific),
  transmute(counts, centre, nW, panel = PANEL_B, series = LAB_TIME,
            value = 100 * timeSpecific / (signals - untested))) %>%
  mutate(panel = factor(panel, c(PANEL_A, PANEL_B)), series = factor(series, c(LAB_ALL, LAB_TIME)))
nullLine <- tibble(panel = factor(PANEL_B, levels(long$panel)), y = 100 * P_EMP)

p <- ggplot(long, aes(centre, value, colour = series)) +
  geom_hline(data = nullLine, aes(yintercept = y), linetype = 'dashed', colour = 'grey40') +
  geom_line(na.rm = TRUE) +
  geom_point(aes(size = nW), na.rm = TRUE) +
  facet_wrap(~ panel, ncol = 1, scales = 'free_y') +
  scale_x_continuous(breaks = seq(0, 21, 3), limits = c(0, 23)) +
  scale_y_continuous(limits = c(0, NA)) +
  scale_colour_manual(values = c('#2b5d8a', '#b2443a'), name = NULL) +
  scale_size_continuous(range = c(1, 3.5), name = 'Window donors') +
  labs(x = 'Window centre (h, DIP time)', y = NULL,
       title = sprintf('%s: %g h rolling windows', TISSUE, WINDOW_WIDTH_H),
       caption = sprintf('Dashed line in B: %g%% expected without a time effect (before the winner\'s curse).', 100 * P_EMP)) +
  theme_bw(base_size = 11) +
  theme(strip.background = element_blank(), strip.text = element_text(hjust = 0, face = 'bold'),
        panel.grid.minor = element_blank(), legend.position = 'bottom', legend.box = 'vertical')

stem <- paste0(FIG_DIR, 'associationsPerWindow', if (P_COL != 'pEmp' || P_EMP != 0.05) sprintf('_%s%g', P_COL, P_EMP) else '')
ggsave(paste0(stem, '.pdf'), p, width = 7, height = 7)
ggsave(paste0(stem, '.png'), p, width = 7, height = 7, dpi = 200)
message('Wrote ', stem, '.{pdf,png} and associationsPerWindow.csv')
