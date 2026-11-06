# =============================================================================
# Shared functions. Sourced by 00_config.R; relies on the paths defined there.
#
#   Data loading     loadDIP, genotypedDonors, gtexDonorId, tissueNames,
#                    selectTissues, loadTissue
#   Windows          dipWindow, permutedRadian, referenceSample
#   Latent factors   runPEER, subsetCovariates
#   Python config    writePythonConfig
# =============================================================================

# --- Data loading --------------------------------------------------------------

# Donor inferred phase (DIP, radians) keyed by GTEx subject ID.
loadDIP <- function() {
  load(DIP_FILE)   # provides `phi`
  data.table(SUBJID = paste0('GTEX-', names(phi)), radian = as.numeric(phi))
}

# Subject IDs of the donors in the genotype .fam file.
genotypedDonors <- function() {
  fam <- fread(paste0(GENOTYPE_PLINK, '.fam'), header = FALSE)
  gtexDonorId(fam$V2)
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

# 24 h reference sample for one window: nW donors drawn at random from all the
# tissue's donors, whatever their DIP, so about 1/N_WINDOWS of them are also
# window donors. Seeded by tissue, permutation and window.
referenceSample <- function(nDonor, nW, perm, tissueIdx, window) {
  set.seed(SEED + 1e6 * tissueIdx + 100 * perm + window)
  sort(sample(nDonor, nW))
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

# Fixed covariates plus PEER recomputed on one donor subset (donors x covariates).
# Covariates constant within the subset are dropped.
subsetCovariates <- function(tis, rows) {
  C <- cbind(tis$cov[rows, , drop = FALSE], runPEER(tis$E[rows, , drop = FALSE]))
  C[, apply(C, 2, \(x) var(x) > 0), drop = FALSE]
}

# --- Python configuration ------------------------------------------------------

# Paths and constants for 03_tensorQTLscan.py, so they are set only in 00_config.R.
writePythonConfig <- function() {
  cfg <- list(GENOTYPE_PLINK = GENOTYPE_PLINK, GTEX_EXPR_DIR = GTEX_EXPR_DIR,
              TISSUE_DIR = TISSUE_DIR, PEER_DIR = PEER_DIR, COUNT_DIR = COUNT_DIR,
              HIT_DIR = HIT_DIR, N_PERM = N_PERM, PERM_BLOCK = PERM_BLOCK,
              N_WINDOWS = N_WINDOWS, CIS_WINDOW = CIS_WINDOW,
              MAF_THRESHOLD = MAF_THRESHOLD, MIN_MAC = MIN_MAC, TQ_NPERM = TQ_NPERM,
              GENE_FDR = GENE_FDR, QVALUE_LAMBDA = QVALUE_LAMBDA,
              VANISH_P = VANISH_P, HET_P = HET_P, SEED = SEED)
  jsonlite::write_json(cfg, paste0(RESULTS_DIR, 'config.json'), auto_unbox = TRUE,
                       digits = NA, pretty = TRUE)
}
