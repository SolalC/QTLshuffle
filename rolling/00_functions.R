# =============================================================================
# Rolling-window functions. Sourced by 00_config.R. Loading, PEER and the
# beta approximation are copied from ../prototype/00_functions.R so that this
# folder runs on its own.
#
#   Data loading   loadDIP, gtexDonorId, tissueDonors, loadTissue, readBim,
#                  readPlinkChromosome, dosage
#   Windows        circularHourDiff, windowRows, referenceSample
#   Covariates     runPEER, windowCovariates, covariateBasis, residualise
#   Association    meanImpute, minorAlleleCount, pFromR2, betaShape1FromDof,
#                  betaApprox, colMaxAbs, geneScan, pairTest
#   Calls          gtexQvalue, clumpGene
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

# Donors of a tissue with expression, covariates, a DIP and a genotype (sorted,
# so that 01 and 02 agree without reading each other's output).
tissueDonors <- function(tissueName) {
  dip <- loadDIP()
  genoDonors <- gtexDonorId(fread(paste0(GENOTYPE_PLINK, '.fam'), header = FALSE)$V2)
  exprIndv <- str_subset(names(fread(exprPath(tissueName), nrows = 0)), '^GTEX')
  covIndv  <- names(fread(covPath(tissueName), nrows = 0))[-1]
  sort(Reduce(intersect, list(dip$SUBJID, exprIndv, covIndv, genoDonors)))
}

# Expression of all genes (donors x genes; PEER uses all of them), gene positions,
# fixed covariates (sex, platform, pcr, PC1-5) and DIP of one tissue.
loadTissue <- function(tissueName) {
  donors <- tissueDonors(tissueName)
  dip  <- loadDIP()
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
       radian = dip$radian[match(donors, dip$SUBJID)] %% (2 * pi),
       E = E, pos = pos, cov = covMat)
}

readBim <- function(prefix) {
  fread(paste0(prefix, '.bim'), header = FALSE,
        col.names = c('chr', 'id', 'cm', 'pos', 'a1', 'a2'))
}

# Genotype codes of the .bim rows `idx` (contiguous, one chromosome) for the
# .fam samples `keep`, read from the SNP-major .bed in chunks. Returns a raw
# matrix (samples x variants): copies of the second .bim allele (0, 1, 2), 3 =
# missing. Raw keeps a whole chromosome in memory at 1 byte per genotype.
readPlinkChromosome <- function(prefix, idx, nSample, keep, chunk = 20000) {
  stopifnot(all(diff(idx) == 1))
  nByte <- ceiling(nSample / 4)
  code  <- c(0L, 3L, 1L, 2L)   # 2-bit codes 00, 01 (missing), 10, 11
  lut   <- sapply(0:255, \(b) code[bitwAnd(bitwShiftR(b, c(0, 2, 4, 6)), 3L) + 1])
  con <- file(paste0(prefix, '.bed'), 'rb'); on.exit(close(con))
  if (!identical(readBin(con, 'raw', 3), as.raw(c(0x6c, 0x1b, 0x01))))
    stop('Not a SNP-major plink .bed')
  G <- matrix(as.raw(3), length(keep), length(idx))
  for (s in seq(1, length(idx), by = chunk)) {
    j <- s:min(s + chunk - 1, length(idx))
    seek(con, 3 + (idx[j[1]] - 1) * nByte)
    b <- as.integer(readBin(con, 'raw', nByte * length(j)))
    codes <- matrix(lut[, b + 1], 4 * nByte, length(j))
    G[, j] <- as.raw(codes[keep, , drop = FALSE])
  }
  G
}

# Raw genotype codes to dosages (NA for missing).
dosage <- function(Graw) {
  G <- matrix(as.numeric(as.integer(Graw)), nrow(Graw), ncol(Graw), dimnames = dimnames(Graw))
  G[G == 3] <- NA
  G
}

# --- Windows -------------------------------------------------------------------

# Signed circular difference (h) between DIP hours and a centre, in [-12, 12).
circularHourDiff <- function(hours, centre) ((hours - centre + 12) %% 24) - 12

# Donors (row indices) whose DIP lies in [centre - W/2, centre + W/2), wrapping
# around midnight (the 0 h window takes 22:30-01:30).
windowRows <- function(radian, centre, width = WINDOW_WIDTH_H) {
  d <- circularHourDiff(radian * 24 / (2 * pi), centre)
  which(d >= -width / 2 & d < width / 2)
}

# 24 h reference sample `ref` for window `centreIdx`: nW donors drawn at random
# from all donors, whatever their DIP (window donors can be drawn).
referenceSample <- function(nDonor, nW, centreIdx, ref) {
  set.seed(SEED + 1e6 + 1e4 * centreIdx + ref)
  sort(sample(nDonor, nW))
}

# --- Covariates ----------------------------------------------------------------

# k PEER factors from the expression of one donor set (donors x all genes).
# 'peer': PEER runs in its own R 3.4 (PEER_RSCRIPT) through peer_fit.R; data are
# exchanged as version-2 .rds files in a temporary folder.
runPEER <- function(E, k) {
  if (PEER_BACKEND == 'peer') {
    tmp <- tempfile('peer_'); dir.create(tmp); on.exit(unlink(tmp, recursive = TRUE))
    inFile <- file.path(tmp, 'Y.rds'); outFile <- file.path(tmp, 'X.rds')
    saveRDS(unname(as.matrix(E)), inFile, version = 2)
    # Clean environment for R 3.4: no R 4.4 R_HOME, library paths or module libraries.
    clean <- c('R_HOME', 'R_LIBS', 'R_LIBS_USER', 'R_LIBS_SITE', 'R_ENVIRON', 'R_ENVIRON_USER',
               'R_PROFILE', 'R_PROFILE_USER', 'R_SHARE_DIR', 'R_INCLUDE_DIR', 'R_DOC_DIR',
               'LD_LIBRARY_PATH')
    status <- system2('env', c(rbind('-u', clean), PEER_RSCRIPT,
                               shQuote(file.path(CODE_DIR, 'peer_fit.R')), shQuote(inFile),
                               shQuote(outFile), k, PEER_ITER, SEED))
    if (status != 0 || !file.exists(outFile)) stop('PEER failed (', PEER_RSCRIPT, ')')
    X <- readRDS(outFile)
  } else if (PEER_BACKEND == 'pca_test_only') {
    X <- prcomp(E, center = TRUE)$x[, seq_len(k), drop = FALSE]
  } else stop('Unknown PEER_BACKEND: ', PEER_BACKEND)
  X <- as.matrix(X)
  dimnames(X) <- list(rownames(E), paste0('PEER', seq_len(ncol(X))))
  X
}

# Fixed covariates plus peerK(n) PEER factors recomputed on the donor set `rows`.
windowCovariates <- function(tis, rows) {
  cbind(tis$cov[rows, , drop = FALSE], runPEER(tis$E[rows, , drop = FALSE], peerK(length(rows))))
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

colMaxAbs <- if (requireNamespace('matrixStats', quietly = TRUE)) {
  \(X) matrixStats::colMaxs(abs(X))
} else {
  \(X) apply(abs(X), 2, max)
}

# cis scan of one gene in one donor set, following tensorQTL map_cis:
# expression and genotypes residualised on the covariates; lead variant = largest
# |correlation|; NPERM_GENO permutations of the residualised expression give the
# null distribution of the gene's minimum p; beta-approximated p with tensorQTL's
# effective degrees of freedom (betaApprox). One deliberate difference: the
# permuted residuals are residualised again (Freedman-Lane), because tensorQTL's
# plain residual permutation is anti-conservative when the covariates are many
# relative to n (~24 covariates for ~90 window donors; see ../prototype/README.md).
# y: expression (n), G: dosages (n x cis variants), Q: covariate basis.
# Returns the gene-level result and every pair with nominal p < PAIR_P_STORE.
geneScan <- function(y, G, Q, seed, nperm = NPERM_GENO) {
  n  <- length(y); df <- n - ncol(Q) - 1
  G  <- meanImpute(G)
  G  <- G[, minorAlleleCount(G) >= MIN_MAC, drop = FALSE]
  if (!ncol(G)) return(NULL)
  Gr <- residualise(G, Q)
  ok <- colSums(Gr^2) > 1e-8
  Gr <- Gr[, ok, drop = FALSE]
  if (!ncol(Gr)) return(NULL)
  yr <- as.vector(residualise(y, Q))
  sxx <- colSums(Gr^2); syy <- sum(yr^2)
  Gn <- sweep(Gr, 2, sqrt(sxx), '/')
  yn <- yr / sqrt(syy)

  r    <- as.vector(crossprod(Gn, yn))
  tval <- r * sqrt(df / (1 - r^2))
  p    <- 2 * pt(-abs(tval), df)
  beta <- r * sqrt(syy / sxx)
  se   <- abs(beta / tval)
  j <- which.max(abs(r))

  set.seed(seed)
  rMax <- unlist(map(split(seq_len(nperm), ceiling(seq_len(nperm) / 1000)), \(chunk) {
    P <- vapply(chunk, \(i) yn[sample.int(n)], numeric(n))
    # Freedman-Lane: permuted residuals residualised on the covariates again.
    colMaxAbs(crossprod(Gn, P) / rep(sqrt(1 - colSums(crossprod(Q, P)^2)), each = ncol(Gn)))
  }))
  ba <- betaApprox(rMax^2, r[j]^2, df)
  keep <- which(p < PAIR_P_STORE)
  list(lead = tibble(variant = colnames(Gr)[j], nVar = ncol(Gr), n = n, dof = df,
                     beta = beta[j], se = se[j], pNominal = p[j],
                     pPerm = (1 + sum(rMax >= abs(r[j]))) / (1 + nperm),
                     pBeta = ba[['pBeta']], betaShape1 = ba[['shape1']],
                     betaShape2 = ba[['shape2']], trueDof = ba[['trueDof']]),
       pairs = tibble(variant = colnames(Gr)[keep], beta = beta[keep], se = se[keep],
                      pNominal = p[keep]))
}

# Variant-gene pairs re-tested in another donor set (same model). G: dosages
# (n x variants, one column per pair); Y: expression (n x genes, one column per
# pair). Pairs with minor allele count < MIN_MAC in the set are NA.
pairTest <- function(Y, G, Q) {
  df <- nrow(Y) - ncol(Q) - 1
  G <- meanImpute(G)
  Gr <- residualise(G, Q); Yr <- residualise(Y, Q)
  sxx <- colSums(Gr^2)
  b  <- colSums(Gr * Yr) / sxx
  se <- sqrt(pmax(colSums(Yr^2) - b^2 * sxx, 0) / df / sxx)
  p  <- 2 * pt(-abs(b / se), df)
  bad <- minorAlleleCount(G) < MIN_MAC | sxx < 1e-8
  b[bad] <- NA; se[bad] <- NA; p[bad] <- NA
  tibble(beta = b, se = se, p = p)
}

# --- Calls -----------------------------------------------------------------------

# eGenes and nominal thresholds as tensorQTL post.calculate_qvalues (GTEx):
# Storey q-value of pBeta (R qvalue, lambda = QVALUE_LAMBDA); eGene if q <=
# GENE_FDR; pthreshold midway between the largest pBeta with q <= GENE_FDR and
# the smallest with q > GENE_FDR; gene-level nominal threshold
# qbeta(pthreshold, shape1, shape2).
gtexQvalue <- function(genes) {
  q <- qvalue::qvalue(genes$pBeta, lambda = QVALUE_LAMBDA)$qvalues
  lb <- max(genes$pBeta[q <= GENE_FDR], -Inf)
  ub <- min(genes$pBeta[q > GENE_FDR], Inf)
  pthreshold <- if (is.finite(lb) && is.finite(ub)) (lb + ub) / 2 else if (is.finite(lb)) lb else NA
  mutate(genes, qval = q, eGene = q <= GENE_FDR, pthreshold = pthreshold,
         pvalNominalThreshold = qbeta(pthreshold, betaShape1, betaShape2))
}

# Greedy clumping of one gene's significant variants (as plink --clump with
# p1 = p2 = the gene's nominal threshold, r2 = CLUMP_R2, no distance limit
# within the cis window): the most significant remaining variant is an index
# variant; variants with r2 > CLUMP_R2 with it join its clump. G: dosages of the
# variants (donors x variants), ordered by increasing p. Returns, per variant,
# the index variant of its clump.
clumpGene <- function(G) {
  G <- meanImpute(G)
  index <- rep(NA_character_, ncol(G)); names(index) <- colnames(G)
  left <- seq_len(ncol(G))
  while (length(left)) {
    i <- left[1]
    r2 <- suppressWarnings(cor(G[, i], G[, left]))^2
    r2[is.na(r2)] <- 0
    member <- left[r2 > CLUMP_R2 | left == i]
    index[member] <- colnames(G)[i]
    left <- setdiff(left, member)
  }
  index
}
