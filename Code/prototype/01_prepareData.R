# =============================================================================
# 01  Prototype data: one tissue (TISSUE), the N_GENES most rhythmic genes
#
#   - donors with expression, covariates, DIP and genotypes;
#   - expression of all genes (PEER is fitted on all of them) and fixed
#     covariates (sex, platform, pcr, PC1-5);
#   - the N_GENES genes with the smallest 24 h harmonic regression p
#     (cd-eQTL step 01);
#   - their cis variants (+/- CIS_WINDOW of the TSS, autosomes, MAF >=
#     MAF_THRESHOLD among the tissue's donors), read from the plink .bed.
#
# Run:     sbatch 01_prepareData.sh
# Usage:   Rscript 01_prepareData.R
# Outputs: Results/prototype/prototypeData.rds, Results/prototype/genes.csv
# =============================================================================

.a <- commandArgs(FALSE); .f <- sub('^--file=', '', .a[grep('^--file=', .a)])
source(file.path(if (length(.f)) dirname(normalizePath(.f)) else '.', '00_config.R'))

fam <- fread(paste0(GENOTYPE_PLINK, '.fam'), header = FALSE)
genoDonors <- gtexDonorId(fam$V2)
tis <- loadTissue(TISSUE, loadDIP(), genoDonors)
message(TISSUE, ': ', length(tis$donors), ' donors, ', ncol(tis$E), ' genes')

stripChr <- \(x) sub('^chr', '', x)
harmonic <- fread(paste0(HARMONIC_DIR, TISSUE, '_LRT_harmonic_linear.csv'))
genes <- tis$pos %>%
  inner_join(harmonic %>% select(gene_id = gene, harmonicP = p.value), by = 'gene_id') %>%
  filter(stripChr(chr) %in% as.character(1:22)) %>%
  arrange(harmonicP, gene_id) %>%
  slice_head(n = N_GENES)

bim <- readBim(GENOTYPE_PLINK)[, row := .I][stripChr(chr) %in% stripChr(genes$chr)]
cisRows <- map(seq_len(nrow(genes)), \(i)
  bim[stripChr(chr) == stripChr(genes$chr[i]) & abs(pos - genes$start[i]) <= CIS_WINDOW, row])
rows <- sort(unique(unlist(cisRows)))
message(length(rows), ' candidate cis variants')

G <- readPlinkVariants(GENOTYPE_PLINK, rows, nrow(fam))[match(tis$donors, genoDonors), , drop = FALSE]
af <- colMeans(G, na.rm = TRUE) / 2
keep <- pmin(af, 1 - af) >= MAF_THRESHOLD
G <- G[, keep, drop = FALSE]; rows <- rows[keep]
variants <- bim[match(rows, row)]
dimnames(G) <- list(tis$donors, variants$id)

genes$nCis <- 0L
cisIdx <- map(seq_len(nrow(genes)), \(i) which(rows %in% cisRows[[i]]))
names(cisIdx) <- genes$gene_id
genes$nCis <- lengths(cisIdx)
genes <- filter(genes, nCis > 0)
cisIdx <- cisIdx[genes$gene_id]
message(nrow(genes), ' genes with cis variants; ', ncol(G), ' variants with MAF >= ', MAF_THRESHOLD)

fwrite(genes, paste0(RESULTS_DIR, 'genes.csv'))
saveRDS(list(tis = tis, genes = genes, G = G, variants = variants, cisIdx = cisIdx),
        DATA_FILE, compress = FALSE)
