# =============================================================================
# 02  Tissue preparation (one tissue per SLURM array task)
#
# Selects the N_TISSUES tissues with the most donors, then for one of them:
#   - expression, fixed covariates (sex, platform, pcr, PC1-5) and DIP;
#   - labels genes rhythmic over 24 h (cd-eQTL step 01, Bonferroni < 0.05);
#   - lists cis variants (+/- CIS_WINDOW of the TSS, tissue MAF >= MAF_THRESHOLD);
#   - computes each gene's effective number of tests (all donors of the tissue).
#
# Run:     sbatch 02_prepareTissues.sh      (array 1-N_TISSUES; after 01)
# Usage:   Rscript 02_prepareTissues.R <tissue index> <n cores>
# Outputs: Results/tissues/selectedTissues.csv, Results/tissues/<tissue>.rds
# =============================================================================

.a <- commandArgs(FALSE); .f <- sub('^--file=', '', .a[grep('^--file=', .a)])
source(file.path(if (length(.f)) dirname(normalizePath(.f)) else '.', '00_config.R'))
library(parallel)

args  <- commandArgs(trailingOnly = TRUE)
i     <- as.numeric(args[1])
nCore <- as.numeric(args[2])

dip <- loadDIP()
genoFiles <- list.files(GENO_DIR, pattern = '^chr[0-9]+\\.rds$')
genoDonors <- rownames(readRDS(paste0(GENO_DIR, genoFiles[1]))$G)
dip <- dip[SUBJID %in% genoDonors]

tissues <- selectTissues(dip)
fwrite(tissues, paste0(TISSUE_DIR, 'selectedTissues.csv'))
tissueName <- tissues$tissue[i]
message('Tissue: ', tissueName, ' (', tissues$n[i], ' donors)')
tis <- loadTissue(tissueName, dip)

# 24 h rhythmicity labels:
harmonic <- fread(paste0(HARMONIC_DIR, tissueName, '_LRT_harmonic_linear.csv'))
rhythmicGenes <- harmonic$gene[p.adjust(harmonic$p.value, 'bonferroni') < 0.05]

# Cis variants and effective number of tests, one chromosome at a time:
tis$genoRow <- match(tis$donors, genoDonors)
tis$cis <- list()
geneInfo <- list()
for (chrom in intersect(unique(tis$pos$chr), str_remove(genoFiles, '\\.rds$'))) {
  geno <- readRDS(paste0(GENO_DIR, chrom, '.rds'))
  G    <- meanImpute(geno$G[tis$genoRow, , drop = FALSE])
  maf  <- colMeans(G) / 2
  okVar <- which(pmin(maf, 1 - maf) >= MAF_THRESHOLD)
  posv <- geno$info$pos[okVar]

  genes <- tis$pos[chr == chrom]
  cisCols <- map2(genes$start, genes$end, \(s, e) {
    lo <- findInterval(s - CIS_WINDOW - 1, posv) + 1
    hi <- findInterval(e + CIS_WINDOW, posv)
    if (hi >= lo) okVar[lo:hi] else integer(0)
  })
  names(cisCols) <- genes$gene_id
  cisCols <- cisCols[lengths(cisCols) > 0]
  if (!length(cisCols)) next

  meff <- unlist(mclapply(cisCols, \(cc) effectiveTests(G[, cc, drop = FALSE]),
                          mc.cores = nCore))
  cols <- sort(unique(unlist(cisCols)))
  tis$cis[[chrom]] <- list(cols = cols, geneIdx = map(cisCols, \(cc) match(cc, cols)))
  geneInfo[[chrom]] <- data.table(gene_id = names(cisCols), chr = chrom,
                                nCis = lengths(cisCols), meff = meff)
  message(chrom, ': ', length(cisCols), ' genes, ', length(cols), ' cis variants')
  rm(geno, G); gc()
}

tis$genes <- rbindlist(geneInfo)[, rhythmic := gene_id %in% rhythmicGenes]
tis$E <- tis$E[, tis$genes$gene_id, drop = FALSE]
message(nrow(tis$genes), ' genes with cis variants; ', sum(tis$genes$rhythmic), ' rhythmic')
saveRDS(tis, paste0(TISSUE_DIR, tissueName, '.rds'), compress = FALSE)
