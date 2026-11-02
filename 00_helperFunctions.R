# =============================================================================
# Shared functions. Sourced by 00_config.R; relies on the paths defined there.
#
#   Data loading     loadDIP, loadGenotypes, gtexDonorId, tissueNames,
#                    selectTissues, loadTissue
#   Windows          dipWindow, permutedRadian
#   Linear algebra   covariateBasis, residualise, meanImpute, minorAlleleCount
#   Multiple tests   effectiveTests
#   Latent factors   runPEER
#   Window scan      cisScan, leadVariantTest, scanWindow
# =============================================================================

# --- Data loading --------------------------------------------------------------

# Donor inferred phase (DIP, radians) keyed by GTEx subject ID.
loadDIP <- function() {
  load(DIP_FILE)   # provides `phi`
  data.table(SUBJID = paste0('GTEX-', names(phi)), radian = as.numeric(phi))
}

# Full plink genotype set (8.25M variants; needs ~100 GB of memory). A .rds path
# holding list(X, bim) in genio's layout is accepted for the toy test.
loadGenotypes <- function() {
  if (grepl('\\.rds$', GENOTYPE_PLINK)) return(readRDS(GENOTYPE_PLINK))
  genio::read_plink(GENOTYPE_PLINK)
}

# WGS sample ID (GTEX-XXXX-....) to subject ID (GTEX-XXXX), as in cd-eQTL.
gtexDonorId <- function(x) paste0('GTEX-', map_chr(str_split(x, '-'), 2))

tissueNames <- function() {
  f <- list.files(GTEX_EXPR_DIR, pattern = '.v10.normalized_expression.bed.gz$')
  str_remove(f, '.v10.normalized_expression.bed.gz')
}

exprPath <- function(t) paste0(GTEX_EXPR_DIR, t, '.v10.normalized_expression.bed.gz')
covPath  <- function(t) paste0(GTEX_COV_DIR, t, '.v10.covariates.txt')

# Donors of a tissue with expression, covariates and a DIP.
tissueDonors <- function(tissueName, dip) {
  exprIndv <- str_subset(names(fread(exprPath(tissueName), nrows = 0)), '^GTEX')
  covIndv  <- names(fread(covPath(tissueName), nrows = 0))[-1]
  Reduce(intersect, list(dip$SUBJID, exprIndv, covIndv))
}

# The N_TISSUES tissues with the most donors (ties broken by name).
selectTissues <- function(dip) {
  tibble(tissue = tissueNames()) %>%
    mutate(n = map_int(tissue, \(t) length(tissueDonors(t, dip)))) %>%
    arrange(desc(n), tissue) %>%
    slice_head(n = N_TISSUES)
}

# Expression (donors x genes), gene positions, fixed covariates (sex, platform,
# pcr, PC1-5; the GTEx PEER factors are not used) and DIP of one tissue.
loadTissue <- function(tissueName, dip) {
  donors <- tissueDonors(tissueName, dip)
  expr <- fread(exprPath(tissueName))
  pos  <- expr[, .(chr = `#chr`, start, end, gene_id)]
  E <- t(as.matrix(expr[, ..donors]))
  dimnames(E) <- list(donors, pos$gene_id)

  cov <- fread(covPath(tissueName))
  covMat <- t(as.matrix(cov[, ..donors]))
  dimnames(covMat) <- list(donors, cov$ID)
  fixed <- intersect(c('sex', 'platform', 'pcr', paste0('PC', 1:5)), colnames(covMat))
  covMat <- covMat[, fixed, drop = FALSE]
  # Drop constant covariates (sex in single-sex tissues):
  covMat <- covMat[, apply(covMat, 2, \(x) length(unique(x)) > 1), drop = FALSE]

  list(name = tissueName, donors = donors,
       radian = dip$radian[match(donors, dip$SUBJID)],
       E = E, pos = pos, cov = covMat)
}

# --- Windows -------------------------------------------------------------------

# Window index (1..N_WINDOWS) of each DIP; window 1 starts at 0 radians.
dipWindow <- function(radian) {
  floor((radian %% (2 * pi)) / (2 * pi / N_WINDOWS)) + 1
}

# DIP for permutation `perm` of tissue `tissueIdx`: perm 0 is the observed DIP;
# otherwise DIP is shuffled across donors, keeping genotype, expression and
# covariates together.
permutedRadian <- function(radian, perm, tissueIdx) {
  if (perm == 0) return(radian)
  set.seed(SEED + 10000 * tissueIdx + perm)
  sample(radian)
}

# --- Linear algebra ------------------------------------------------------------

# Orthonormal basis of [1, covariates], dropping columns that are constant or
# collinear within the subset.
covariateBasis <- function(C) {
  C <- C[, apply(C, 2, \(x) var(x) > 0), drop = FALSE]
  q <- qr(cbind(1, C))
  qr.Q(q)[, seq_len(q$rank), drop = FALSE]
}

residualise <- function(Y, Q) Y - Q %*% crossprod(Q, Y)

meanImpute <- function(G) {
  storage.mode(G) <- 'double'
  na <- which(is.na(G), arr.ind = TRUE)
  if (nrow(na)) G[na] <- colMeans(G, na.rm = TRUE)[na[, 2]]
  G
}

minorAlleleCount <- function(G) {
  ac <- colSums(G)
  pmin(ac, 2 * nrow(G) - ac)
}

# --- Multiple testing ----------------------------------------------------------

# Effective number of independent tests among a gene's cis variants (eigenMT-
# style): within consecutive blocks of MEFF_BLOCK variants, the number of
# eigenvalues of the correlation matrix needed to explain MEFF_VAR of the
# variance, summed over blocks. Computed once per gene on all donors.
effectiveTests <- function(G) {
  blocks <- split(seq_len(ncol(G)), ceiling(seq_len(ncol(G)) / MEFF_BLOCK))
  sum(map_dbl(blocks, \(b) {
    g <- G[, b, drop = FALSE]
    g <- g[, apply(g, 2, var) > 0, drop = FALSE]
    if (ncol(g) <= 1) return(ncol(g))
    ev <- pmax(eigen(cor(g), symmetric = TRUE, only.values = TRUE)$values, 0)
    which(cumsum(ev) / sum(ev) >= MEFF_VAR)[1]
  }))
}

# --- Latent factors ------------------------------------------------------------

# N_PEER factors from the expression of one donor subset (donors x genes).
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

# --- Window scan ---------------------------------------------------------------

# Covariate basis of a donor subset: fixed covariates plus PEER recomputed on the
# subset's expression.
subsetBasis <- function(tis, rows) {
  covariateBasis(cbind(tis$cov[rows, , drop = FALSE], runPEER(tis$E[rows, , drop = FALSE])))
}

# Marginal regression of every cis variant on residualised expression, one
# chromosome at a time. Returns the lead variant per gene (lowest p):
# gene, chr, variant column (in the chromosome file), beta, se, p, nTested.
cisScan <- function(tis, geno, rows, Q) {
  n  <- length(rows)
  df <- n - ncol(Q) - 1
  Er <- residualise(tis$E[rows, , drop = FALSE], Q)
  map_dfr(names(tis$cis), \(chr) {
    cc   <- tis$cis[[chr]]
    Graw <- meanImpute(geno[[chr]]$G[tis$genoRow[rows], cc$cols, drop = FALSE])
    ok   <- minorAlleleCount(Graw) >= MIN_MAC
    Gr   <- residualise(Graw[, ok, drop = FALSE], Q)
    sxx  <- colSums(Gr^2)
    okMap <- rep(NA_integer_, length(ok)); okMap[ok] <- seq_len(sum(ok))
    map_dfr(names(cc$geneIdx), \(gene) {
      idx  <- cc$geneIdx[[gene]]          # positions in cc$cols
      m    <- okMap[idx]                  # positions in Gr
      keep <- !is.na(m)
      keep[keep] <- sxx[m[keep]] > 1e-8
      if (!any(keep)) return(NULL)
      idx <- idx[keep]; m <- m[keep]
      e   <- Er[, gene]
      b   <- drop(crossprod(Gr[, m, drop = FALSE], e)) / sxx[m]
      rss <- pmax(sum(e^2) - b^2 * sxx[m], 0)
      se  <- sqrt(rss / df / sxx[m])
      p   <- 2 * pt(-abs(b / se), df)
      j   <- which.min(p)
      tibble(gene = gene, chr = chr, col = cc$cols[idx[j]],
             beta = b[j], se = se[j], p = p[j], nTested = length(m))
    })
  })
}

# Test the lead variant of each hit in another donor subset (same model).
leadVariantTest <- function(tis, geno, rows, Q, hits) {
  n  <- length(rows)
  df <- n - ncol(Q) - 1
  Er <- residualise(tis$E[rows, hits$gene, drop = FALSE], Q)
  pmap_dfr(list(hits$gene, hits$chr, hits$col), \(gene, chr, col) {
    g <- meanImpute(geno[[chr]]$G[tis$genoRow[rows], col, drop = FALSE])
    if (minorAlleleCount(g) < MIN_MAC)
      return(tibble(beta_out = NA_real_, se_out = NA_real_, p_out = NA_real_))
    gr  <- residualise(g, Q)
    sxx <- sum(gr^2)
    e   <- Er[, gene]
    b   <- sum(gr * e) / sxx
    se  <- sqrt(max(sum(e^2) - b^2 * sxx, 0) / df / sxx)
    tibble(beta_out = b, se_out = se, p_out = 2 * pt(-abs(b / se), df))
  })
}

# One window: eGenes among the window donors (Meff-corrected lead p, BH across
# genes), then the lead variant of each eGene re-tested on the out-of-window
# donors. A vanishing hit has out-of-window p > VANISH_P and a significant
# window vs out-of-window difference (z test, p < HET_P).
# Returns list(counts = one row, hits = hit-level table).
scanWindow <- function(tis, geno, inWindow) {
  rowsW <- which(inWindow); rowsO <- which(!inWindow)
  base  <- tibble(nW = length(rowsW), nO = length(rowsO))
  empty <- base %>% mutate(nGenes = NA_integer_, nHits = NA_integer_,
                           nVanish = NA_integer_, nHitsRhythmic = NA_integer_,
                           nVanishRhythmic = NA_integer_, nOutUntestable = NA_integer_)
  if (length(rowsW) < MIN_WINDOW_N) return(list(counts = empty, hits = NULL))

  lead <- cisScan(tis, geno, rowsW, subsetBasis(tis, rowsW)) %>%
    left_join(tis$genes %>% select(gene = gene_id, meff, rhythmic), by = 'gene') %>%
    mutate(pGene = pmin(1, p * meff), qGene = p.adjust(pGene, 'BH'))
  hits <- lead %>% filter(qGene < GENE_FDR)

  if (nrow(hits)) {
    hits <- bind_cols(hits, leadVariantTest(tis, geno, rowsO, subsetBasis(tis, rowsO), hits)) %>%
      mutate(z = (beta - beta_out) / sqrt(se^2 + se_out^2),
             pHet = 2 * pnorm(-abs(z)),
             vanish = !is.na(p_out) & p_out > VANISH_P & pHet < HET_P)
  } else {
    hits <- hits %>% mutate(beta_out = numeric(), se_out = numeric(), p_out = numeric(),
                            z = numeric(), pHet = numeric(), vanish = logical())
  }

  counts <- base %>% mutate(
    nGenes          = nrow(lead),
    nHits           = nrow(hits),
    nVanish         = sum(hits$vanish),
    nHitsRhythmic   = sum(hits$rhythmic),
    nVanishRhythmic = sum(hits$vanish & hits$rhythmic),
    nOutUntestable  = sum(is.na(hits$p_out)))
  list(counts = counts, hits = hits)
}
