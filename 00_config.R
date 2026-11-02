# =============================================================================
# Paths and constants shared by every R script. Each script sources this file
# first; it in turn sources 00_helperFunctions.R.
#
# Every input path can be overridden with an environment variable of the same
# name prefixed by QTLSHUFFLE_ (used by tests/run_toy_test.sh). Defaults point
# to the Bunya locations used by the cd-eQTL analysis.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
})

.envDir <- function(var, default) {
  d <- Sys.getenv(var, default)
  if (!endsWith(d, '/')) d <- paste0(d, '/')
  d
}
.envPath <- function(var, default) Sys.getenv(var, default)
.envNum  <- function(var, default) as.numeric(Sys.getenv(var, default))

ROOT        <- .envDir('QTLSHUFFLE_ROOT', '/scratch/project_mnt/S0007/solal/QTLshuffle/')
RESULTS_DIR <- .envDir('QTLSHUFFLE_RESULTS', paste0(ROOT, 'Results/'))

# --- Inputs (same files as the cd-eQTL analysis) ---------------------------------
CIRCADIAN_ROOT  <- '/scratch/project_mnt/S0007/solal/circadian/'
# GTEx v10 normalised expression (inverse-normal transformed across all donors)
GTEX_EXPR_DIR   <- .envDir('QTLSHUFFLE_GTEX_EXPR_DIR', '/QRISdata/Q8106/Public/QTL/eQTL_expression_matrices/')
# GTEx v10 covariates; only sex, platform, pcr and PC1-5 are used (PEER is recomputed)
GTEX_COV_DIR    <- .envDir('QTLSHUFFLE_GTEX_COV_DIR', '/QRISdata/Q8106/Public/QTL/eQTL_covariates/')
# GTEx WGS genotypes (plink). A .rds path is accepted for the toy test only.
GENOTYPE_PLINK  <- .envPath('QTLSHUFFLE_GENOTYPE_PLINK', paste0(
  '/QRISdata/Q8106/Controlled/Genotype/filtered/',
  'GTEx_Analysis_2021-02-11_v9_WholeGenomeSeq_953Indiv_cQTL_subset_inferred'))
# Donor inferred phase (CHIRAL, Talamanca et al. 2023): R object `phi`, radians
DIP_FILE        <- .envPath('QTLSHUFFLE_DIP_FILE', paste0(CIRCADIAN_ROOT, 'Data/CHIRAL/DIP.RData'))
# 24 h rhythmicity of every gene (cd-eQTL step 01), used only to label genes
HARMONIC_DIR    <- .envDir('QTLSHUFFLE_HARMONIC_DIR', paste0(CIRCADIAN_ROOT, 'Results/inferred/harmonicReg/'))
# Genes used by CHIRAL to infer DIP, one Ensembl ID per line (no header).
# TODO: provide this file; while it is missing no genes are excluded.
CHIRAL_GENES    <- .envPath('QTLSHUFFLE_CHIRAL_GENES', paste0(ROOT, 'Data/CHIRAL_genes.txt'))

# --- Intermediate results ------------------------------------------------------
GENO_DIR    <- paste0(RESULTS_DIR, 'genotypes/')     # 01: per-chromosome dosages
TISSUE_DIR  <- paste0(RESULTS_DIR, 'tissues/')       # 02: tissue objects
COUNT_DIR   <- paste0(RESULTS_DIR, 'windowCounts/')  # 03: counts per window and permutation
HIT_DIR     <- paste0(RESULTS_DIR, 'windowHits/')    # 03: hit-level results, observed data only
SUMMARY_DIR <- paste0(RESULTS_DIR, 'summary/')       # 04
for (.d in c(GENO_DIR, TISSUE_DIR, COUNT_DIR, HIT_DIR, SUMMARY_DIR))
  dir.create(.d, recursive = TRUE, showWarnings = FALSE)

# --- Analysis constants --------------------------------------------------------
N_TISSUES     <- .envNum('QTLSHUFFLE_N_TISSUES', 10)  # tissues with the most donors
N_WINDOWS     <- 8          # non-overlapping 3 h DIP windows (2*pi / 8 radians each)
MIN_WINDOW_N  <- 60         # windows with fewer donors are skipped (NA counts)
CIS_WINDOW    <- 1e6        # bp either side of the TSS
MAF_THRESHOLD <- 0.01       # minimum MAF across the tissue's donors
MIN_MAC       <- 10         # minimum minor allele count within a window / out-of-window set
N_PEER        <- 10         # PEER factors, recomputed within every donor subset
PEER_ITER     <- 1000       # PEER maximum iterations
MEFF_BLOCK    <- 200        # eigenMT-style blocks of consecutive cis variants
MEFF_VAR      <- 0.995      # share of LD-matrix variance retained per block
GENE_FDR      <- 0.05       # BH FDR across genes on Meff-corrected lead p-values
VANISH_P      <- 0.05       # out-of-window lead-variant p above this ...
HET_P         <- 0.05       # ... and window vs out-of-window heterogeneity p below this
N_PERM        <- .envNum('QTLSHUFFLE_N_PERM', 1000)   # permutations of DIP (0 = observed)
PERM_BLOCK    <- .envNum('QTLSHUFFLE_PERM_BLOCK', 20) # permutations per SLURM task
SEED          <- 20261102
SEX_SPECIFIC  <- c('Ovary', 'Prostate', 'Testis', 'Uterus', 'Vagina')

# PEER backend. 'peer' (the R package) is the analysis. 'pca_test_only' exists
# solely so tests/run_toy_test.sh can run where PEER is not installed.
PEER_BACKEND  <- .envPath('QTLSHUFFLE_PEER_BACKEND', 'peer')

# --- Shared functions ----------------------------------------------------------
.a <- commandArgs(FALSE); .f <- sub('^--file=', '', .a[grep('^--file=', .a)])
CODE_DIR <- if (length(.f)) dirname(normalizePath(.f)) else '.'
source(file.path(CODE_DIR, '00_helperFunctions.R'))
