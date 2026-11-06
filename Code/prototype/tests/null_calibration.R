# Null calibration of geneScan: 2000 simulated genes with no genetic effect
# (n = 80, 19 covariates, real toy genotypes), Freedman-Lane permutations vs
# the plain residual permutations used by tensorQTL. Run after run_toy_test.sh.
Sys.setenv(QTLSHUFFLE_RESULTS='tests/toy/Results/', QTLSHUFFLE_NPERM_GENO='1000', QTLSHUFFLE_PEER_BACKEND='pca_test_only')
source('00_config.R'); library(parallel)
d <- readRDS(DATA_FILE)
src <- paste(readLines('00_functions.R'), collapse = '\n')
plain <- sub('rep(sqrt(1 - colSums(crossprod(Q, P)^2)), each = ncol(Gn))', '1', src, fixed = TRUE)
stopifnot(plain != src)
e <- new.env(); eval(parse(text = plain), envir = e); geneScanPlain <- e$geneScan
set.seed(7); rows <- sample(nrow(d$G), 80)
Q <- covariateBasis(cbind(d$tis$cov[rows, ], matrix(rnorm(80 * 10), 80)))
run <- function(f) unlist(mclapply(1:2000, \(i) {
  set.seed(1e6 + i); y <- rnorm(80)
  f(y, d$G[rows, d$cisIdx[[(i %% 50) + 1]], drop = FALSE], Q, i)$pBeta }, mc.cores = 4))
for (nm in c('freedmanLane', 'plainResidual')) {
  p <- run(if (nm == 'freedmanLane') geneScan else geneScanPlain)
  r <- sapply(c(0.05, 0.01, 0.001), \(a) mean(p < a))
  cat(sprintf('%s: n genes %d; P(p<0.05)=%.4f [SE %.4f]; P(p<0.01)=%.4f [SE %.4f]; P(p<0.001)=%.4f; KS p=%.3g\n',
      nm, length(p), r[1], sqrt(.05*.95/2000), r[2], sqrt(.01*.99/2000), r[3], ks.test(p, 'punif')$p.value))
}
