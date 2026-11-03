# =============================================================================
# 04  eQTL effect against inferred phase (observed DIP; no permutations)
#
# For each selected gene, one cis variant is followed across the day. Two
# choices of variant, one figure each:
#   24hLead  the lead variant of a standard 24 h scan (all donors, largest
#            marginal effect): what a standard eQTL analysis reports;
#   GxTLead  the variant with the strongest genotype x time interaction
#            (cosinor F test, all donors): the candidate cd-eQTL. A variant
#            with a pure phase or amplitude effect has no 24 h marginal effect,
#            so the 24 h scan does not pick it.
# For the chosen variant:
#   - effect (beta, SE) re-estimated in each of the N_WINDOWS non-overlapping
#     3 h DIP windows, PEER recomputed on the window donors (as in 02); one
#     point per window, at the window centre;
#   - cosinor G x time fit on all donors,
#     y ~ g + g:cos(DIP) + g:sin(DIP) + cos(DIP) + sin(DIP) + covariates,
#     drawn as the expected window effect b + lambda (g1 cos t + g2 sin t),
#     lambda = sinc(w / 2) for a window of w radians (0.97 for 3 h).
# The variant is chosen on all donors, which include the window donors, so the
# points are not independent of the choice. This figure is descriptive.
#
# Run:     sbatch 04_phasePlot.sh            (after 01; 1 + N_WINDOWS PEER fits)
# Outputs: Results_prototype/phasePlot/{windowBeta,cosinorFit}.csv,
#          betaByPhase_24hLead.pdf, betaByPhase_GxTLead.pdf
# =============================================================================

.a <- commandArgs(FALSE); .f <- sub('^--file=', '', .a[grep('^--file=', .a)])
source(file.path(if (length(.f)) dirname(normalizePath(.f)) else '.', '00_config.R'))

outDir <- paste0(RESULTS_DIR, 'phasePlot/')
dir.create(outDir, showWarnings = FALSE)
d <- readRDS(DATA_FILE)
radian  <- d$tis$radian %% (2 * pi)
cosr    <- cos(radian); sinr <- sin(radian)
allRows <- seq_along(radian)
genes   <- d$genes$gene_id
cisG    <- \(gene, rows) d$G[rows, d$cisIdx[[gene]], drop = FALSE]

# --- Variant choice on all donors ----------------------------------------------
CAll  <- subsetCovariates(d$tis, allRows)
QAll  <- covariateBasis(CAll)
QTime <- covariateBasis(cbind(CAll, cosr, sinr))
leads <- bind_rows(
  map_dfr(genes, \(gene) {
    res <- geneScan(d$tis$E[, gene], cisG(gene, allRows), QAll, seed = 0, nperm = 0)
    if (!is.null(res)) tibble(gene = gene, variant = res$variant, nVar = res$nVar)
  }) %>% mutate(choice = '24hLead'),
  map_dfr(genes, \(gene) {
    res <- gxtLead(d$tis$E[, gene], cisG(gene, allRows), QTime, cosr, sinr)
    tibble(gene = gene, variant = res$variant, nVar = res$nVar)
  }) %>% mutate(choice = 'GxTLead'))

# --- Cosinor G x time fit on all donors -------------------------------------------
width  <- 2 * pi / N_WINDOWS
lambda <- sin(width / 2) / (width / 2)
cosinor <- leads %>% mutate(fit = map2(gene, variant, \(gene, variant) {
  y <- d$tis$E[, gene]
  g <- as.vector(meanImpute(matrix(d$G[, variant])))
  full <- lm(y ~ g * (cosr + sinr) + CAll)
  cf <- coef(full)
  tibble(beta24h = coef(summary(lm(y ~ g + CAll)))['g', 1],
         b = cf[['g']], g1 = cf[['g:cosr']], g2 = cf[['g:sinr']],
         pGxT = anova(lm(y ~ g + cosr + sinr + CAll), full)$`Pr(>F)`[2])
})) %>% unnest(fit)

# --- Effect in each window (PEER per window) ----------------------------------------
win <- dipWindow(radian)
windowBeta <- map_dfr(seq_len(N_WINDOWS), \(w) {
  rows <- which(win == w)
  if (length(rows) < MIN_WINDOW_N) return(NULL)
  Q <- covariateBasis(subsetCovariates(d$tis, rows))
  map_dfr(seq_len(nrow(leads)), \(i) {
    x <- leadVariantTest(d$tis$E[rows, leads$gene[i]], d$G[rows, leads$variant[i]], Q)
    tibble(choice = leads$choice[i], gene = leads$gene[i], variant = leads$variant[i],
           window = w, centre = (w - 0.5) * width, n = length(rows),
           beta = x[['beta']], se = x[['se']], p = x[['p']])
  })
})

fwrite(cosinor, paste0(outDir, 'cosinorFit.csv'))
fwrite(windowBeta, paste0(outDir, 'windowBeta.csv'))

# --- Figures: one panel per gene, ordered by the G x time p ---------------------------
toHours <- 24 / (2 * pi)
for (ch in c('24hLead', 'GxTLead')) {
  cz  <- filter(cosinor, choice == ch) %>% arrange(pGxT) %>%
    mutate(panel = factor(sprintf('%s  %s\nGxT p = %.1e', gene, variant, pGxT)) %>% fct_inorder())
  key <- select(cz, gene, panel)
  curve <- cz %>% crossing(t = seq(0, 2 * pi, length.out = 100)) %>%
    mutate(beta = b + lambda * (g1 * cos(t) + g2 * sin(t)))
  pts <- filter(windowBeta, choice == ch) %>% left_join(key, by = 'gene')
  p <- ggplot(pts, aes(centre * toHours, beta)) +
    geom_hline(yintercept = 0, colour = 'grey85') +
    geom_hline(data = cz, aes(yintercept = beta24h), linetype = 'dashed', colour = 'grey50') +
    geom_line(data = curve, aes(t * toHours, beta), colour = 'steelblue') +
    geom_linerange(aes(ymin = beta - 1.96 * se, ymax = beta + 1.96 * se), linewidth = 0.3, colour = 'grey40') +
    geom_linerange(aes(ymin = beta - se, ymax = beta + se), linewidth = 1) +
    geom_point(size = 1) +
    facet_wrap(~panel, scales = 'free_y', ncol = 5) +
    scale_x_continuous(breaks = seq(0, 24, 6), limits = c(0, 24)) +
    labs(title = paste0(TISSUE, ': ', ifelse(ch == '24hLead', '24 h lead variant', 'strongest G x time variant')),
         x = 'Inferred phase (DIP, h)', y = 'Variant effect (beta)',
         caption = paste0('Points: one 3 h DIP window each, PEER recomputed per window; ',
                          'thick bar = +/- 1 SE (SD of the estimate), thin bar = 95% CI. ',
                          'Dashed: 24 h effect. Blue: cosinor G x time fit on all donors (window-attenuated).')) +
    theme_classic(base_size = 8) +
    theme(strip.background = element_blank())
  ggsave(paste0(outDir, 'betaByPhase_', ch, '.pdf'), p, width = 10,
         height = 2 * ceiling(nrow(cz) / 5) + 0.7, limitsize = FALSE)
}
message('Wrote ', outDir)
