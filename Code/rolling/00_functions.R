# =============================================================================
# Rolling-window functions. Sourced by 00_config.R. Loading and PEER are
# copied from ../prototype/00_functions.R so that this folder runs on its own. Genotypes are read with the BEDMatrix and genio
# packages.
#
#   Data loading   loadDIP, gtexDonorId, genotypeSampleIds, tissueDonors, loadTissue,
#                  openGenotypes, readGenotypes (BEDMatrix, genio)
#   Windows        circularHourDiff, windowRows, referenceSample
#   Covariates     runPEER, windowCovariates
#   Association    meanImpute, minorAlleleCount, pairTest (lm)
# The window scan uses MatrixEQTL (03) and clumping uses plink2 --clump (04).
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

# Sample IDs of the genotypes: the plink .fam once 02a has made it (genio),
# otherwise the #CHROM line of the VCF header (grep stops at the first match).
genotypeSampleIds <- function() {
  if (file.exists(paste0(GENOTYPE_PLINK, '.fam')))
    return(genio::read_fam(GENOTYPE_PLINK, verbose = FALSE)$id)
  header <- fread(cmd = paste('zcat', shQuote(GENOTYPE_VCF), "| grep -m 1 '^#CHROM'"),
                  sep = '\t', header = TRUE, nrows = 0)
  names(header)[-(1:9)]
}

# Donors of a tissue with expression, covariates, a DIP and a genotype (sorted,
# so that the scripts agree without reading each other's output).
tissueDonors <- function(tissueName) {
  dip <- loadDIP()
  genoDonors <- gtexDonorId(genotypeSampleIds())
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

# Genotypes of the plink file GENOTYPE_PLINK (made by 02a), memory-mapped with
# BEDMatrix: nothing is read until a subset is taken. `rows` maps `donors` (GTEx
# subject IDs) to .fam rows; `bim` has chr (without 'chr'), id, pos, alt (A1)
# and ref (A2). BEDMatrix counts copies of the A1 allele (alt for a plink2
# conversion); see tests/check_bedmatrix.R.
openGenotypes <- function(donors) {
  fam <- genio::read_fam(GENOTYPE_PLINK, verbose = FALSE)
  bim <- genio::read_bim(GENOTYPE_PLINK, verbose = FALSE)
  bim$chr <- sub('^chr', '', bim$chr)
  if (anyDuplicated(bim$id)) stop('Duplicated variant IDs in ', GENOTYPE_PLINK, '.bim')
  rows <- match(donors, gtexDonorId(fam$id))
  if (anyNA(rows)) stop(sum(is.na(rows)), ' donors missing from ', GENOTYPE_PLINK, '.fam')
  list(X = BEDMatrix::BEDMatrix(GENOTYPE_PLINK, n = nrow(fam), p = nrow(bim)),
       rows = rows, bim = bim)
}

# Dosages (NA if missing) of donors `donorIdx` (indices into `donors`) at .bim
# rows `cols`, as a numeric matrix with variant IDs as column names.
readGenotypes <- function(geno, donorIdx, cols) {
  G <- geno$X[geno$rows[donorIdx], cols, drop = FALSE]
  G <- matrix(as.numeric(G), nrow(G), ncol(G))
  colnames(G) <- geno$bim$id[cols]
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

# --- Association -----------------------------------------------------------------

meanImpute <- function(G) {
  na <- which(is.na(G), arr.ind = TRUE)
  if (nrow(na)) G[na] <- colMeans(G, na.rm = TRUE)[na[, 2]]
  G
}

minorAlleleCount <- function(G) {
  ac <- colSums(G)
  pmin(ac, 2 * nrow(G) - ac)
}

# Variant-gene pairs re-tested in another donor set with lm():
# y ~ g + covariates, one model per pair (columns j of Y and G). C: covariate
# matrix of the donor set (fixed + PEER). Missing genotypes are mean-imputed,
# as in the scan; pairs with minor allele count < MIN_MAC in the set are NA.
pairTest <- function(Y, G, C) {
  map_dfr(seq_len(ncol(Y)), \(j) {
    g <- G[, j]
    g[is.na(g)] <- mean(g, na.rm = TRUE)
    cf <- if (min(sum(g), 2 * length(g) - sum(g)) >= MIN_MAC) coef(summary(lm(Y[, j] ~ g + C)))
    if (is.null(cf) || !'g' %in% rownames(cf)) return(tibble(beta = NA_real_, se = NA_real_, p = NA_real_))
    tibble(beta = cf['g', 1], se = cf['g', 2], p = cf['g', 4])
  })
}


# --- Annotation (figures) ----------------------------------------------------------------

# Gene symbol (GENCODE gene_name) for versioned or unversioned Ensembl IDs;
# the Ensembl ID is kept where the GTF has no entry.
geneSymbol <- function(genes, gtf = GENE_GTF) {
  map <- read_tsv(gtf, comment = '#', col_names = FALSE, col_types = 'ccciicccc', progress = FALSE) %>%
    filter(X3 == 'gene') %>%
    transmute(id = sub('\\..*', '', str_match(X9, 'gene_id "([^"]+)"')[, 2]),
              symbol = str_match(X9, 'gene_name "([^"]+)"')[, 2]) %>%
    distinct(id, .keep_all = TRUE)
  s <- map$symbol[match(sub('\\..*', '', genes), map$id)]
  coalesce(s, genes)
}

# rsID for pipeline variant IDs (chr:pos:REF:ALT, 02a) from the GTEx lookup
# table (variant_id chr_pos_REF_ALT_b38); the pipeline ID is kept where the
# table has no rsID. Only the requested rows are read (grep).
variantRsid <- function(variants, lookup = VARIANT_LOOKUP) {
  gtexId <- paste0(gsub(':', '_', variants), '_b38')
  ids <- tempfile(); writeLines(c('variant_id', unique(gtexId)), ids)
  tab <- fread(cmd = sprintf('zcat %s | grep -wF -f %s', shQuote(lookup), ids),
               select = c('variant_id', 'rs_id_dbSNP155_GRCh38p13'))
  unlink(ids)
  rs <- tab$rs_id_dbSNP155_GRCh38p13[match(gtexId, tab$variant_id)]
  if_else(is.na(rs) | rs == '.', variants, rs)
}
