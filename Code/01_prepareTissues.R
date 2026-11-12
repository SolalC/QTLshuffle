# =============================================================================
# 01  Tissue preparation (one tissue per SLURM array task)
#
# Selects the N_TISSUES tissues with the most donors that have a genotype and a
# DIP, then for one of them saves expression, fixed covariates (sex, platform,
# pcr, PC1-5) and DIP for the PEER step, and labels every gene as rhythmic over
# 24 h or not (cd-eQTL step 01, Bonferroni < 0.05). Also writes the
# configuration read by 03_tensorQTLscan.py.
#
# Run:     sbatch 01_prepareTissues.sh      (array 1-N_TISSUES)
# Usage:   Rscript 01_prepareTissues.R <tissue index>
# Outputs: Results/config.json, Results/main/tissues/selectedTissues.csv,
#          Results/tissues/<tissue>.rds, Results/main/tissues/<tissue>_genes.csv
# =============================================================================

.a <- commandArgs(FALSE); .f <- sub('^--file=', '', .a[grep('^--file=', .a)])
source(file.path(if (length(.f)) dirname(normalizePath(.f)) else '.', '00_config.R'))

i <- as.numeric(commandArgs(trailingOnly = TRUE)[1])

writePythonConfig()
dip <- loadDIP()
dip <- dip[SUBJID %in% genotypedDonors()]
message(nrow(dip), ' donors with genotypes and a DIP')

tissues <- selectTissues(dip)
fwrite(tissues, paste0(TISSUE_DIR, 'selectedTissues.csv'))
tissueName <- tissues$tissue[i]
message('Tissue: ', tissueName, ' (', tissues$n[i], ' donors)')
tis <- loadTissue(tissueName, dip)

harmonic <- fread(paste0(HARMONIC_DIR, tissueName, '_LRT_harmonic_linear.csv'))
rhythmicGenes <- harmonic$gene[p.adjust(harmonic$p.value, 'bonferroni') < 0.05]
genes <- data.table(gene_id = tis$pos$gene_id, rhythmic = tis$pos$gene_id %in% rhythmicGenes)
message(nrow(genes), ' genes; ', sum(genes$rhythmic), ' rhythmic')

fwrite(genes, paste0(TISSUE_DIR, tissueName, '_genes.csv'))
saveRDS(tis, paste0(TISSUE_DIR, tissueName, '.rds'), compress = FALSE)
