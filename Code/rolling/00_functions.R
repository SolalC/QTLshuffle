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
#   Scan and call  scanChromosome (MatrixEQTL), callSignals (BH + plink2 --clump),
#                  pairBeta (MatrixEQTL beta of given pairs)
# The same scan and call functions are used for windows (03, 04) and for
# reference samples treated as windows (07), so that both are selected alike.
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


# --- Scan and call (windows 03/04, pseudo-windows 07, re-tests 08) ------------------

slicedData <- function(M) { s <- SlicedData$new(); s$CreateFromMatrix(M); s$ResliceCombined(2000); s }

# cis scan of one donor set x one chromosome with MatrixEQTL (modelLINEAR):
#   y ~ g + covariates, every variant within CIS_WINDOW of the TSS (`end` of the
#   expression .bed) with minor allele count >= MIN_MAC in the donor set
#   (missing genotypes mean-imputed). rows: donor indices into tis$donors;
#   cov: donors x covariates. Returns pairs with p < pThreshold (gene,
#   variant, beta, se, t, pNominal), genes (gene, nVar, minP), the number of
#   cis tests and of variants tested.
scanChromosome <- function(tis, geno, rows, cov, chrom, pThreshold = PAIR_P_STORE) {
  cov   <- cov[, apply(cov, 2, var) > 0, drop = FALSE]     # e.g. one platform in the set
  genes <- tis$pos[sub('^chr', '', chr) == as.character(chrom)]
  cols  <- which(geno$bim$chr == as.character(chrom))
  G <- meanImpute(readGenotypes(geno, rows, cols))
  ok <- minorAlleleCount(G) >= MIN_MAC
  G <- G[, ok, drop = FALSE]; pos <- geno$bim$pos[cols][ok]
  me <- Matrix_eQTL_main(
    snps = slicedData(t(G)), gene = slicedData(t(tis$E[rows, genes$gene_id, drop = FALSE])),
    cvrt = slicedData(t(cov)), useModel = modelLINEAR, errorCovariance = numeric(),
    output_file_name = tempfile(), pvOutputThreshold = 0,             # no trans tests
    output_file_name.cis = tempfile(), pvOutputThreshold.cis = pThreshold,
    snpspos = data.frame(snp = colnames(G), chr = chrom, pos = pos),
    genepos = data.frame(gene = genes$gene_id, chr = chrom, left = genes$end, right = genes$end),
    cisDist = CIS_WINDOW, pvalue.hist = FALSE, min.pv.by.genesnp = TRUE,
    noFDRsaveMemory = FALSE, verbose = FALSE)
  # cis variants per gene (as MatrixEQTL defines cis: |pos - TSS| <= cisDist)
  nVar <- findInterval(genes$end + CIS_WINDOW, pos) - findInterval(genes$end - CIS_WINDOW - 1, pos)
  list(pairs = as_tibble(me$cis$eqtls) %>%
         transmute(gene = as.character(gene), variant = as.character(snps),
                   beta, se = beta / statistic, t = statistic, pNominal = pvalue),
       genes = tibble(gene = genes$gene_id, nVar = nVar, minP = me$cis$min.pv.gene[genes$gene_id]) %>%
         filter(nVar > 0),
       nTests = me$cis$ntests, nVariants = ncol(G), nCov = ncol(cov))
}

# plink2 --clump of one set's significant variants of one gene; returns the
# index variants and clump sizes. LD from all the tissue's donors (PLINK_CHR).
clumpGene <- function(x) {
  dir <- tempfile('clump_'); dir.create(dir); on.exit(unlink(dir, recursive = TRUE))
  fwrite(tibble(ID = x$variant, P = pmax(x$pNominal, 1e-300)), file.path(dir, 'assoc.tsv'), sep = '\t')
  status <- system2(PLINK2, c('--bfile', PLINK_CHR(x$chr[1]), '--clump', file.path(dir, 'assoc.tsv'),
                              '--clump-p1', 1, '--clump-p2', 1, '--clump-r2', CLUMP_R2,
                              '--clump-kb', CLUMP_KB, '--threads', 1, '--memory', 4000,
                              '--out', file.path(dir, 'out')), stdout = FALSE, stderr = FALSE)
  out <- file.path(dir, 'out.clumps')
  if (status != 0 || !file.exists(out)) stop('plink2 --clump failed for ', x$set[1], ', ', x$gene[1])
  fread(out) %>% transmute(variant = ID, clumpSize = TOTAL + 1L)
}

# Significant pairs and independent signals of one or more donor sets:
#   BH q-value of every stored pair over all cis tests of its set (nTests:
#   set, nTestsSet), significant pair q < PAIR_FDR, then each significant
#   gene's variants clumped with plink2 (clumpGene). pairs: set, gene, chr,
#   variant, beta, se, t, pNominal. Returns list(sig, independent).
callSignals <- function(pairs, nTests, nCore) {
  pairs <- pairs %>% left_join(nTests, by = 'set') %>%
    group_by(set) %>% mutate(qval = p.adjust(pNominal, 'BH', n = nTestsSet[1])) %>% ungroup()
  sig <- filter(pairs, qval < PAIR_FDR) %>% arrange(set, gene, pNominal)
  if (!nrow(sig)) return(list(sig = sig, independent = mutate(sig, clumpSize = integer(), signal = integer(), nSignal = integer())))
  maxP <- max(sig$pNominal)
  if (maxP > PAIR_P_STORE / 10)
    warning('Significant pairs reach p = ', signif(maxP, 2), ', close to PAIR_P_STORE (', PAIR_P_STORE, ')')
  groups <- split(sig, interaction(sig$set, sig$gene, drop = TRUE))
  index <- parallel::mclapply(groups, \(x) mutate(clumpGene(x), set = x$set[1], gene = x$gene[1]),
                              mc.cores = nCore, mc.preschedule = FALSE)
  failed <- map_lgl(index, \(x) inherits(x, 'try-error'))
  if (any(failed)) stop(sum(failed), ' clumping runs failed, e.g. ', as.character(index[[which(failed)[1]]]))
  independent <- bind_rows(index) %>%
    inner_join(sig, by = c('set', 'gene', 'variant')) %>%
    arrange(set, gene, pNominal) %>%
    group_by(set, gene) %>% mutate(signal = row_number(), nSignal = n()) %>% ungroup()
  list(sig = sig, independent = independent)
}

# MatrixEQTL beta and SE of given variant-gene pairs (gene, variant; one
# chromosome) in one donor set, same model and filters as scanChromosome.
# Pairs whose variant has minor allele count < MIN_MAC in the set are NA.
pairBeta <- function(tis, geno, rows, cov, pairs, chrom) {
  cov   <- cov[, apply(cov, 2, var) > 0, drop = FALSE]
  genes <- tis$pos[match(unique(pairs$gene), gene_id)]
  cols  <- match(unique(pairs$variant), geno$bim$id)
  G <- meanImpute(readGenotypes(geno, rows, cols))
  ok <- minorAlleleCount(G) >= MIN_MAC
  out <- pairs %>% select(gene, variant) %>% mutate(beta = NA_real_, se = NA_real_)
  if (!any(ok)) return(out)
  G <- G[, ok, drop = FALSE]; pos <- geno$bim$pos[cols][ok]
  me <- Matrix_eQTL_main(
    snps = slicedData(t(G)), gene = slicedData(t(tis$E[rows, genes$gene_id, drop = FALSE])),
    cvrt = slicedData(t(cov)), useModel = modelLINEAR, errorCovariance = numeric(),
    output_file_name = tempfile(), pvOutputThreshold = 0,
    output_file_name.cis = tempfile(), pvOutputThreshold.cis = 1,   # every cis pair
    snpspos = data.frame(snp = colnames(G), chr = chrom, pos = pos),
    genepos = data.frame(gene = genes$gene_id, chr = chrom, left = genes$end, right = genes$end),
    cisDist = CIS_WINDOW, pvalue.hist = FALSE, noFDRsaveMemory = FALSE, verbose = FALSE)
  res <- as_tibble(me$cis$eqtls) %>%
    transmute(gene = as.character(gene), variant = as.character(snps), b = beta, s = beta / statistic)
  out %>% left_join(res, by = c('gene', 'variant')) %>% transmute(gene, variant, beta = b, se = s)
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
