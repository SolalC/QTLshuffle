# =============================================================================
# 03  cis-eQTL scan of one window centre x one chromosome (MatrixEQTL)
#
# Window donors only. For every gene on the chromosome (TSS = `end` of the
# expression .bed), every cis variant within CIS_WINDOW of the TSS with minor
# allele count >= MIN_MAC among the window donors (missing genotypes
# mean-imputed) is tested with MatrixEQTL's linear model:
#   y ~ g + sex + platform + pcr + PC1-5 + peerK(n) PEER factors
# (window covariates from 01). Nominal p-values only: no permutations. Pairs
# with p < PAIR_P_STORE are stored, with the number of cis tests, for the
# FDR in 04. Genotypes are read from the plink file of 02a (BEDMatrix).
#
# Run:     sbatch 03_windowScan.sh       (array 1..24 x length(CHROMOSOMES))
# Usage:   Rscript 03_windowScan.R <task>
#          task = (centre index - 1) * length(CHROMOSOMES) + chromosome index
# Outputs: Results/rolling/w<W>h/scan/c<centre>_chr<c>_{pairs,genes}.csv.gz
#            pairs: gene, variant, beta, se, t, pNominal (p < PAIR_P_STORE)
#            genes: gene, nVar (cis variants tested), minP; nTests in every row
# =============================================================================

.a <- commandArgs(FALSE); .f <- sub('^--file=', '', .a[grep('^--file=', .a)])
source(file.path(if (length(.f)) dirname(normalizePath(.f)) else '.', '00_config.R'))
task <- as.integer(commandArgs(trailingOnly = TRUE)[1])
centreIdx <- (task - 1) %/% length(CHROMOSOMES) + 1
chrom     <- CHROMOSOMES[(task - 1) %% length(CHROMOSOMES) + 1]
centre    <- CENTRES_H[centreIdx]
outFile   <- \(what) sprintf('%sc%02d_chr%d_%s.csv.gz', SCAN_DIR, centre, chrom, what)

win <- readRDS(WINDOW_FILE)
if (!win$windows$used[centreIdx]) {
  message('Centre ', centre, ' h has ', win$windows$n[centreIdx], ' donors (< ', MIN_WINDOW_N, '): skipped')
  fwrite(tibble(), outFile('genes')); fwrite(tibble(), outFile('pairs'))
  quit(save = 'no')
}
tis  <- readRDS(TISSUE_FILE)
rows <- win$windows$rows[[centreIdx]]
cov  <- win$cov[[centreIdx]]
cov  <- cov[, apply(cov, 2, var) > 0, drop = FALSE]   # e.g. one platform in the window

genes <- tis$pos[sub('^chr', '', chr) == as.character(chrom)]
geno  <- openGenotypes(tis$donors)
cols  <- which(geno$bim$chr == as.character(chrom))
G <- meanImpute(readGenotypes(geno, rows, cols))       # window donors x chromosome variants
rownames(G) <- tis$donors[rows]
ok <- minorAlleleCount(G) >= MIN_MAC
G <- G[, ok, drop = FALSE]; pos <- geno$bim$pos[cols][ok]
message('Centre ', centre, ' h, chr', chrom, ': ', length(rows), ' donors, ', ncol(cov),
        ' covariates, ', nrow(genes), ' genes, ', ncol(G), ' variants with MAC >= ', MIN_MAC)

sliced <- \(M) { s <- SlicedData$new(); s$CreateFromMatrix(M); s$ResliceCombined(2000); s }
t0 <- Sys.time()
me <- Matrix_eQTL_main(
  snps = sliced(t(G)), gene = sliced(t(tis$E[rows, genes$gene_id, drop = FALSE])),
  cvrt = sliced(t(cov)), useModel = modelLINEAR, errorCovariance = numeric(),
  output_file_name = tempfile(), pvOutputThreshold = 0,             # no trans tests
  output_file_name.cis = tempfile(), pvOutputThreshold.cis = PAIR_P_STORE,
  snpspos = data.frame(snp = colnames(G), chr = chrom, pos = pos),
  genepos = data.frame(gene = genes$gene_id, chr = chrom, left = genes$end, right = genes$end),
  cisDist = CIS_WINDOW, pvalue.hist = FALSE, min.pv.by.genesnp = TRUE,
  noFDRsaveMemory = FALSE, verbose = FALSE)
message('MatrixEQTL: ', me$cis$ntests, ' cis tests in ',
        round(as.numeric(difftime(Sys.time(), t0, units = 'mins')), 1), ' min')

# cis variants per gene (as MatrixEQTL defines cis: |pos - TSS| <= cisDist)
nVar <- findInterval(genes$end + CIS_WINDOW, pos) - findInterval(genes$end - CIS_WINDOW - 1, pos)
fwrite(tibble(centre = centre, chr = chrom, gene = genes$gene_id, nVar = nVar,
              minP = me$cis$min.pv.gene[genes$gene_id], nTests = me$cis$ntests) %>%
         filter(nVar > 0), outFile('genes'))
fwrite(as_tibble(me$cis$eqtls) %>%
         transmute(centre = centre, gene = as.character(gene), variant = as.character(snps),
                   beta, se = beta / statistic, t = statistic, pNominal = pvalue),
       outFile('pairs'))
