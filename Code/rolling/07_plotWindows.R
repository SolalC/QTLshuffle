# =============================================================================
# 07  Associations per window centre
#
# Panel A: number of independent signals (04, BH q < PAIR_FDR, clumped) per
#          window centre.
# Panel B: number of those signals whose window beta lies outside the
#          distribution of the N_REF random 24 h reference samples
#          (06: pEmp < P_EMP, two-sided). The dashed line is the count expected
#          if no signal depended on time (P_EMP x signals). Selection on the
#          window beta (winner's curse) inflates panel B above that line even
#          without any time effect, so the line is a floor, not a null.
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
         expected = P_EMP * (signals - untested)) %>%
  arrange(centre)
fwrite(counts, paste0(FIG_DIR, 'associationsPerWindow.csv'))
print(as_tibble(counts), n = Inf)

long <- bind_rows(
  transmute(counts, centre, nW, panel = 'A  Independent signals (window q < 0.05)', count = signals, expected = NA_real_),
  transmute(counts, centre, nW, panel = sprintf('B  Signals outside 24 h reference (%s < %g)', P_COL, P_EMP),
            count = timeSpecific, expected))

p <- ggplot(long, aes(centre, count)) +
  geom_line(colour = 'grey60') +
  geom_point(aes(size = nW), colour = '#2b5d8a') +
  geom_line(aes(y = expected), linetype = 'dashed', colour = '#b2443a', na.rm = TRUE) +
  facet_wrap(~ panel, ncol = 1, scales = 'free_y') +
  scale_x_continuous(breaks = seq(0, 21, 3), limits = c(0, 23)) +
  scale_y_continuous(limits = c(0, NA)) +
  scale_size_continuous(range = c(1, 3.5), name = 'Window donors') +
  labs(x = 'Window centre (h, DIP time)', y = 'Number of associations',
       title = sprintf('%s: %g h rolling windows', TISSUE, WINDOW_WIDTH_H),
       caption = sprintf('Dashed line in B: %g x tested signals (expected without time effect, before winner\'s curse).', P_EMP)) +
  theme_bw(base_size = 11) +
  theme(strip.background = element_blank(), strip.text = element_text(hjust = 0, face = 'bold'),
        panel.grid.minor = element_blank())

stem <- paste0(FIG_DIR, 'associationsPerWindow', if (P_COL != 'pEmp' || P_EMP != 0.05) sprintf('_%s%g', P_COL, P_EMP) else '')
ggsave(paste0(stem, '.pdf'), p, width = 7, height = 6)
ggsave(paste0(stem, '.png'), p, width = 7, height = 6, dpi = 200)
message('Wrote ', stem, '.{pdf,png} and associationsPerWindow.csv')
