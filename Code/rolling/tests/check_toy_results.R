# Checks on the rolling-window toy run (toy data from ../tests/make_toy_data.R).
suppressPackageStartupMessages({ library(tidyverse); library(data.table) })
res <- commandArgs(trailingOnly = TRUE)[1]
fail <- function(...) { message('FAIL: ', ...); quit(status = 1) }
w <- file.path(res, 'w3h')

windows <- fread(file.path(w, 'windows.csv'))
if (nrow(windows) != 24) fail('expected 24 window centres')
print(windows)
if (file.exists(file.path(w, 'results', 'missingReference.csv'))) fail('missing reference samples')
out <- fread(file.path(w, 'results', 'rollingWindow.csv.gz'))
if (!nrow(out)) fail('no independent signals')
if (sum(startsWith(names(out), 'beta_ref_')) != 40) fail('expected 40 beta_ref columns')
if (any(out$nRef < 30)) fail('signals with < 30 non-missing reference betas')

# Window donor set check: the 0 h window wraps around midnight.
tis <- readRDS(file.path(res, 'data', 'tissue.rds'))
h <- tis$radian * 12 / pi
if (sum(h >= 22.5 | h < 1.5) != windows[centre == 0, n]) fail('0 h window does not wrap around midnight')

# Toy genes 1-10 per chromosome are phase QTL, 11-20 amplitude QTL, 21-30 mean
# eQTL, 31-40 rhythmic only, 41-150 null. Each has one causal variant, so one
# independent signal per eGene is expected (more = clumping too lenient).
arch <- \(gene) { k <- as.integer(str_sub(gene, 7, 9))
  case_when(k <= 10 ~ 'phase', k <= 20 ~ 'amplitude', k <= 30 ~ 'mean', k <= 40 ~ 'rhythmic only', TRUE ~ 'null') }
summ <- out %>% mutate(arch = arch(gene)) %>% group_by(arch) %>%
  summarise(signals = n(), centreGene = n_distinct(paste(centre, gene)),
            meanSignalsPerEGene = signals / centreGene,
            pEmpBelow05 = mean(pEmp < 0.05), pEmpLargerBelow05 = mean(pEmpLarger < 0.05))
print(summ)
mean_ <- filter(summ, arch == 'mean'); tv <- filter(summ, arch %in% c('phase', 'amplitude'))
if (sum(tv$signals) == 0) fail('no phase or amplitude signals')
if (weighted.mean(tv$pEmpLargerBelow05, tv$signals) <= 2 * mean_$pEmpLargerBelow05)
  fail('time-varying eQTL not more often outside the 24 h distribution than mean eQTL')
nullRate <- sum(out$gene %in% out$gene[arch(out$gene) %in% c('null', 'rhythmic only')]) / nrow(out)
message('Share of signals in null or rhythmic-only genes: ', round(nullRate, 3))

# Selection-bias control (07, 08): every centre calibrated; calibration must not
# call mean (time-invariant) eQTL time dependent more often than pEmp < 0.05 does.
calSummary <- rbindlist(map(list.files(file.path(w, 'null', 'calibration'), pattern = '_summary\\.csv$', full.names = TRUE), fread))
if (nrow(calSummary) != 24) fail('expected 24 08 summaries, found ', nrow(calSummary))
print(calSummary[, .(centre, nW, windowSignals, pseudoSignalsMean, pctWindowPEmp05, pctNullMean, pWindow, timeDependentCal)])
cal <- list.files(file.path(w, 'null', 'calibration'), pattern = '_window\\.csv\\.gz$', full.names = TRUE)
if (!length(cal)) fail('no centre calibrated')
cal <- rbindlist(map(cal, fread))
if (!setequal(cal$centre, calSummary[windowSignals > 0 & pseudoSignalsMean > 0, centre])) fail('calibrated centres do not match the summaries')
nul <- rbindlist(map(list.files(file.path(w, 'null', 'calibration'), pattern = '_null\\.csv\\.gz$', full.names = TRUE), fread))
if (max(abs(nul$betaSelf - nul$beta), na.rm = TRUE) > 1e-6) fail('pairBeta does not reproduce the pseudo-window scan beta')
calSumm <- cal %>% mutate(arch = arch(gene)) %>% group_by(arch) %>%
  summarise(signals = n(), pEmpKBelow05 = mean(pEmpK < 0.05, na.rm = TRUE),
            pCalBelow05 = mean(pCal < 0.05, na.rm = TRUE), timeDependent = mean(timeDependent, na.rm = TRUE))
print(calSumm)
message('Pseudo-window signals: ', nrow(nul), '; % pEmp < 0.05 among them: ', round(100 * mean(nul$pEmp < 0.05, na.rm = TRUE), 1))
m <- filter(calSumm, arch == 'mean')
if (nrow(m) && m$pCalBelow05 > m$pEmpKBelow05 + 0.02) fail('calibration raises the time-dependent rate of mean eQTL')
message('PASS')
