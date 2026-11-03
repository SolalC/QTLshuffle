# =============================================================================
# Prototype functions. Sourced by 00_config.R.
#
#   Data loading   loadDIP, gtexDonorId, tissueDonors, loadTissue,
#                  readBim, readPlinkVariants
#   Windows        dipWindow, permutedRadian, referenceSample
#   Covariates     runPEER, subsetCovariates, covariateBasis, residualise
#   Association    meanImpute, minorAlleleCount, pFromR2, betaShape1FromDof,
#                  betaApprox, geneScan, gxtLead, leadVariantTest
#   Window scan    geneQvalue, scanWindow
# =============================================================================

# --- Data loading --------------------------------------------------------------

# Donor inferred phase (DIP, radians) keyed by GTEx subject ID.
loadDIP <- function() {
  load(DIP_FILE)   # provides `phi`
  data.table(SUBJID = paste0('GTEX-', names(phi)), radian = as.numeric(phi))
}

# WGS sample ID (GTEX-XXXX-....) to subject ID (GTEX-XXXX), as in cd-eQTL.
gtexDonorId <- function(x) paste0('GTEX-', map_chr(str_split(x, '-'), 2))

exprPath <- function(t) paste0(GTEX_EXPR_DIR, t, '.v10.normalized_expression.bed.gz')
covPath  <- function(t) paste0(GTEX_COV_DIR, t, '.v10.covariates.txt')

# Donors of a tissue with expression, covariates, a DIP and a genotype.
tissueDonors <- function(tissueName, dip, genoDonors) {
  exprIndv <- str_subset(names(fread(exprPath(tissueName), nrows = 0)), '^GTEX')
  covIndv  <- names(fread(covPath(tissueName), nrows = 0))[-1]
  Reduce(intersect, list(dip$SUBJID, exprIndv, covIndv, genoDonors))
}

# Expression of all genes (donors x genes; PEER uses all of them), gene positions,
# fixed covariates (sex, platform, pcr, PC1-5) and DIP of one tissue.
loadTissue <- function(tissueName, dip, genoDonors) {
  donors <- tissueDonors(tissueName, dip, genoDonors)
  expr <- fread(exprPath(tissueName))
  pos  <- expr[, .(chr = `#chr`, start, end, gene_id)]
  E <- t(as.matrix(expr[, ..donors]))
  dimnames(E) <- list(donors, pos$gene_id)

  cov <- fread(covPath(tissueName))
  covMat <- t(as.matrix(cov[, ..donors]))
  dimnames(covMat) <- list(donors, cov$ID)
  fixed <- intersect(c('sex', 'platform', 'pcr', paste0('PC', 1:5)), colnames(covMat))
  covMat <- covMat[, fixed, drop = FALSE]
  covMat <- covMat[, apply(covMat, 2, \(x) length(unique(x)) > 1), drop = FALSE]

  list(name = tissueName, donors = donors,
       radian = dip$radian[match(donors, dip$SUBJID)],
       E = E, pos = pos, cov = covMat)
}

readBim <- function(prefix) {
  fread(paste0(prefix, '.bim'), header = FALSE,
        col.names = c('chr', 'id', 'cm', 'pos', 'a1', 'a2'))
}

# Dosage (copies of the second bim allele; NA if missing) of the variants at
# rows `idx` of the .bim, for all samples of the .fam (samples x variants).
# Reads only those variants from the SNP-major .bed.
readPlinkVariants <- function(prefix, idx, nSample) {
  nByte <- ceiling(nSample / 4)
  code  <- c(0, NA, 1, 2)   # 2-bit codes 00, 01, 10, 11
  lut   <- sapply(0:255, \(b) code[bitwAnd(bitwShiftR(b, c(0, 2, 4, 6)), 3L) + 1])
  con <- file(paste0(prefix, '.bed'), 'rb'); on.exit(close(con))
  magic <- readBin(con, 'raw', 3)
  if (!identical(magic, as.raw(c(0x6c, 0x1b, 0x01)))) stop('Not a SNP-major plink .bed')
  G <- matrix(NA_real_, nSample, length(idx))
  for (j in seq_along(idx)) {
    seek(con, 3 + (idx[j] - 1) * nByte)
    b <- as.integer(readBin(con, 'raw', nByte))
    G[, j] <- as.vector(lut[, b + 1])[seq_len(nSample)]
  }
  G
}

# --- Windows -------------------------------------------------------------------

# Window index (1..N_WINDOWS) of each DIP; window 1 starts at 0 radians.
dipWindow <- function(radian) floor((radian %% (2 * pi)) / (2 * pi / N_WINDOWS)) + 1

# DIP for permutation `perm`: perm 0 is the observed DIP; otherwise DIP is
# shuffled across donors, keeping genotype, expression and covariates together.
permutedRadian <- function(radian, perm) {
  if (perm == 0) return(radian)
  set.seed(SEED + perm)
  sample(radian)
}

# 24 h reference sample for one window: nW donors drawn at random from all
# donors whatever their DIP (about 1/N_WINDOWS are also window donors).
referenceSample <- function(nDonor, nW, perm, window) {
  set.seed(SEED + 1e6 + 100 * perm + window)
  sort(sample(nDonor, nW))
}

# --- Covariates ----------------------------------------------------------------

# N_PEER factors from the expression of one donor set (donors x all genes).
runPEER <- function(E, k = N_PEER) {
  if (PEER_BACKEND == 'peer') {
    model <- peer::PEER()
    peer::PEER_setPhenoMean(model, as.matrix(E))
    peer::PEER_setNk(model, k)
    peer::PEER_setNmax_iterations(model, PEER_ITER)
    peer::PEER_update(model)
    X <- peer::PEER_getX(model)
  } else if (PEER_BACKEND == 'pca_test_only') {
    X <- prcomp(E, center = TRUE)$x[, seq_len(k), drop = FALSE]
  } else stop('Unknown PEER_BACKEND: ', PEER_BACKEND)
  colnames(X) <- paste0('PEER', seq_len(ncol(X)))
  X
}

# Fixed covariates plus PEER recomputed on the donor set `rows`.
subsetCovariates <- function(tis, rows) {
  cbind(tis$cov[rows, , drop = FALSE], runPEER(tis$E[rows, , drop = FALSE]))
}

# Orthonormal basis of [1, covariates], dropping constant or collinear columns.
covariateBasis <- function(C) {
  C <- C[, apply(C, 2, \(x) var(x) > 0), drop = FALSE]
  q <- qr(cbind(1, C))
  qr.Q(q)[, seq_len(q$rank), drop = FALSE]
}

residualise <- function(Y, Q) Y - Q %*% crossprod(Q, Y)

# --- Association ---------------------------------------------------------------

meanImpute <- function(G) {
  na <- which(is.na(G), arr.ind = TRUE)
  if (nrow(na)) G[na] <- colMeans(G, na.rm = TRUE)[na[, 2]]
  G
}

minorAlleleCount <- function(G) {
  ac <- colSums(G)
  pmin(ac, 2 * nrow(G) - ac)
}

# Beta approximation of the permutation p-value, as tensorQTL
# (core.fit_beta_parameters / calculate_beta_approx_pval):
#   1. effective degrees of freedom: the dof at which the permutation minimum
#      p-values, converted from the maximum r2, have a moment-matched first
#      beta shape parameter of 1 (root found on log dof);
#   2. beta shape parameters by maximum likelihood on those p-values (start at
#      the moment estimates);
#   3. p_beta = beta CDF of the nominal p recomputed with the effective dof.
pFromR2 <- function(r2, dof) 2 * pt(-sqrt(dof * r2 / (1 - r2)), dof)

betaShape1FromDof <- function(r2, dof) {
  p <- pFromR2(r2, dof); m <- mean(p); v <- mean((p - m)^2)
  m * (m * (1 - m) / v - 1)
}

betaApprox <- function(r2Perm, r2Nominal, dofInit) {
  f <- \(lx) log(betaShape1FromDof(r2Perm, exp(lx)))
  trueDof <- tryCatch(exp(uniroot(f, log(dofInit) + c(-3, 3), extendInt = 'yes', tol = 1e-4)$root),
                      error = \(e) optimize(\(x) abs(betaShape1FromDof(r2Perm, x) - 1),
                                             c(1, 10 * dofInit))$minimum)
  p <- pFromR2(r2Perm, trueDof); m <- mean(p); v <- mean((p - m)^2)
  s1 <- m * (m * (1 - m) / v - 1); s2 <- s1 * (1 / m - 1)
  nll <- \(s) if (any(s <= 0)) Inf else -sum(dbeta(p, s[1], s[2], log = TRUE))
  sh <- optim(c(s1, s2), nll)$par
  pTrue <- pFromR2(r2Nominal, trueDof)
  c(pBeta = pbeta(pTrue, sh[1], sh[2]), shape1 = sh[1], shape2 = sh[2], trueDof = trueDof)
}

# cis scan of one gene in one donor set, following the tensorQTL map_cis logic:
# expression and genotypes residualised on the covariates; lead variant = largest
# |correlation|; NPERM_GENO permutations of the residualised expression give the
# null distribution of the gene's minimum p; beta-approximated p from a beta
# fitted to it, with tensorQTL's effective degrees of freedom (betaApprox).
# One deliberate difference: tensorQTL correlates the permuted residuals as
# they are, which are no longer orthogonal to the covariates. Their correlations
# then have variance ~1/(n-1) instead of 1/(n-k-1), so the null is too narrow
# and p-values are too small when k covariates are many relative to n (here
# ~19 covariates for ~100 window donors). The permuted residuals are therefore
# residualised again (Freedman-Lane), which restores the null variance.
# y: expression (n), G: dosages (n x cis variants), Q: covariate basis;
# nperm = 0 returns the nominal lead variant only.
geneScan <- function(y, G, Q, seed, nperm = NPERM_GENO) {
  n  <- length(y); df <- n - ncol(Q) - 1
  G  <- meanImpute(G)
  G  <- G[, minorAlleleCount(G) >= MIN_MAC, drop = FALSE]
  if (!ncol(G)) return(NULL)
  Gr <- residualise(G, Q)
  ok <- colSums(Gr^2) > 1e-8
  G <- G[, ok, drop = FALSE]; Gr <- Gr[, ok, drop = FALSE]
  if (!ncol(Gr)) return(NULL)
  yr <- as.vector(residualise(y, Q))
  Gn <- sweep(Gr, 2, sqrt(colSums(Gr^2)), '/')
  yn <- yr / sqrt(sum(yr^2))

  r <- as.vector(crossprod(Gn, yn))
  j <- which.max(abs(r))
  rToP <- \(r) 2 * pt(-abs(r) * sqrt(df / (1 - r^2)), df)
  sxx  <- sum(Gr[, j]^2)
  beta <- sum(Gr[, j] * yr) / sxx
  se   <- sqrt(max(sum(yr^2) - beta^2 * sxx, 0) / df / sxx)

  if (nperm == 0)   # nominal lead variant only
    return(tibble(variant = colnames(G)[j], nVar = ncol(G), beta = beta, se = se,
                  pNominal = rToP(r[j])))
  set.seed(seed)
  rMax <- unlist(map(split(seq_len(nperm), ceiling(seq_len(nperm) / 1000)), \(chunk) {
    P <- vapply(chunk, \(i) yn[sample.int(n)], numeric(n))
    # Freedman-Lane: the permuted residuals are residualised on the covariates
    # again before computing correlations (see note above).
    rP <- crossprod(Gn, P) / rep(sqrt(1 - colSums(crossprod(Q, P)^2)), each = ncol(Gn))
    apply(abs(rP), 2, max)
  }))
  ba <- betaApprox(rMax^2, r[j]^2, df)
  tibble(variant = colnames(G)[j], nVar = ncol(G), beta = beta, se = se,
         pNominal = rToP(r[j]),
         pPerm = (1 + sum(rMax >= abs(r[j]))) / (1 + nperm),
         pBeta = ba[['pBeta']], betaShape1 = ba[['shape1']], betaShape2 = ba[['shape2']],
         trueDof = ba[['trueDof']])
}

# Variant with the strongest genotype x time interaction (cosinor) among a gene's
# cis variants: F test of g:cos + g:sin given g, cos, sin and covariates.
# Q0: covariate basis that already includes cos(DIP) and sin(DIP).
gxtLead <- function(y, G, Q0, cosr, sinr) {
  G  <- meanImpute(G)
  G  <- G[, minorAlleleCount(G) >= MIN_MAC, drop = FALSE]
  yr <- as.vector(residualise(y, Q0))
  df <- length(y) - ncol(Q0) - 3
  f <- vapply(seq_len(ncol(G)), \(j) {
    X <- residualise(cbind(G[, j], G[, j] * cosr, G[, j] * sinr), Q0)
    rss0 <- sum(lm.fit(X[, 1, drop = FALSE], yr)$residuals^2)
    rss1 <- sum(lm.fit(X, yr)$residuals^2)
    ((rss0 - rss1) / 2) / (rss1 / df)
  }, numeric(1))
  j <- which.max(f)
  tibble(variant = colnames(G)[j], nVar = ncol(G), pGxTLead = pf(f[j], 2, df, lower.tail = FALSE))
}

# Lead variant re-tested in another donor set (same model).
leadVariantTest <- function(y, g, Q) {
  df <- length(y) - ncol(Q) - 1
  g <- meanImpute(matrix(g))
  if (minorAlleleCount(g) < MIN_MAC) return(c(beta = NA, se = NA, p = NA))
  gr <- residualise(g, Q); yr <- residualise(y, Q)
  sxx <- sum(gr^2)
  b  <- sum(gr * yr) / sxx
  se <- sqrt(max(sum(yr^2) - b^2 * sxx, 0) / df / sxx)
  c(beta = b, se = se, p = 2 * pt(-abs(b / se), df))
}

# --- Window scan ---------------------------------------------------------------

# Storey q-values (R qvalue, lambda as GTEx). With few genes the pi0 estimate
# can be 0 (no p >= lambda); qvalue then fails and BH (pi0 = 1) is used instead.
geneQvalue <- function(p) {
  q <- tryCatch(qvalue::qvalue(p, lambda = QVALUE_LAMBDA)$qvalues, error = \(e) NULL)
  if (is.null(q)) return(list(q = p.adjust(p, 'BH'), method = 'BH'))
  list(q = q, method = 'storey')
}

# One window: eGenes among the window donors, then each eGene's lead variant
# re-tested in the 24 h reference sample (rowsR). Vanishing hit: reference
# p > VANISH_P and window vs reference difference z test p < HET_P.
scanWindow <- function(d, rowsW, rowsR, seedBase) {
  base <- tibble(nW = length(rowsW), n24h = length(rowsR),
                 nShared = length(intersect(rowsW, rowsR)))
  QW <- covariateBasis(subsetCovariates(d$tis, rowsW))
  lead <- map_dfr(seq_len(nrow(d$genes)), \(i) {
    gene <- d$genes$gene_id[i]
    res <- geneScan(d$tis$E[rowsW, gene], d$G[rowsW, d$cisIdx[[gene]], drop = FALSE], QW,
                    seedBase + i)
    if (!is.null(res)) mutate(res, gene = gene, .before = 1)
  })
  qv <- geneQvalue(lead$pBeta)
  lead <- mutate(lead, qval = qv$q, qMethod = qv$method)
  hits <- filter(lead, qval < GENE_FDR)

  QR <- covariateBasis(subsetCovariates(d$tis, rowsR))
  retest <- map_dfr(seq_len(nrow(hits)), \(h) {
    x <- leadVariantTest(d$tis$E[rowsR, hits$gene[h]], d$G[rowsR, hits$variant[h]], QR)
    tibble(beta_24h = x[['beta']], se_24h = x[['se']], p_24h = x[['p']])
  })
  if (!nrow(hits)) retest <- tibble(beta_24h = numeric(), se_24h = numeric(), p_24h = numeric())
  hits <- bind_cols(hits, retest) %>%
    mutate(z = (beta - beta_24h) / sqrt(se^2 + se_24h^2),
           pHet = 2 * pnorm(-abs(z)),
           vanish = !is.na(p_24h) & p_24h > VANISH_P & pHet < HET_P)

  counts <- base %>% mutate(nGenes = nrow(lead), nHits = nrow(hits),
                            nVanish = sum(hits$vanish),
                            nUntestable24h = sum(is.na(hits$p_24h)),
                            qMethod = qv$method)
  list(counts = counts, hits = hits)
}
