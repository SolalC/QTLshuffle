# =============================================================================
# Prototype (R only): one tissue, the N_GENES most rhythmic genes.
# Paths and constants for every prototype script; sources 00_functions.R.
# Every input path can be overridden by an environment variable prefixed
# QTLSHUFFLE_ (used by tests/run_prototype_toy_test.sh).
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

ROOT        <- .envDir('QTLSHUFFLE_ROOT', '/scratch/user/uqschauq/QTLshuffle/')
RESULTS_DIR <- .envDir('QTLSHUFFLE_RESULTS', paste0(ROOT, 'Results/prototype/'))

# --- Inputs (same files as the cd-eQTL analysis) ---------------------------------
CIRCADIAN_ROOT <- '/scratch/project_mnt/S0007/solal/circadian/'
# GTEx v10 normalised expression (inverse-normal transformed across all donors)
GTEX_EXPR_DIR  <- .envDir('QTLSHUFFLE_GTEX_EXPR_DIR', '/QRISdata/Q8106/Public/QTL/eQTL_expression_matrices/')
# GTEx v10 covariates; only sex, platform, pcr and PC1-5 are used (PEER is recomputed)
GTEX_COV_DIR   <- .envDir('QTLSHUFFLE_GTEX_COV_DIR', '/QRISdata/Q8106/Public/QTL/eQTL_covariates/')
# GTEx WGS genotypes (plink prefix: .bed/.bim/.fam)
GENOTYPE_PLINK <- .envPath('QTLSHUFFLE_GENOTYPE_PLINK', paste0(
  '/QRISdata/Q8106/Controlled/Genotype/filtered/',
  'GTEx_Analysis_2021-02-11_v9_WholeGenomeSeq_953Indiv_cQTL_subset_inferred'))
# Donor inferred phase (CHIRAL): R object `phi`, radians
DIP_FILE       <- .envPath('QTLSHUFFLE_DIP_FILE', paste0(CIRCADIAN_ROOT, 'Data/CHIRAL/DIP.RData'))
# 24 h rhythmicity of every gene (cd-eQTL step 01): ranks genes for selection
HARMONIC_DIR   <- .envDir('QTLSHUFFLE_HARMONIC_DIR', paste0(CIRCADIAN_ROOT, 'Results/inferred/harmonicReg/'))

# --- Results ---------------------------------------------------------------------
DATA_FILE   <- paste0(RESULTS_DIR, 'prototypeData.rds')  # 01
COUNT_DIR   <- paste0(RESULTS_DIR, 'windowCounts/')      # 02: counts per window and DIP permutation
HIT_DIR     <- paste0(RESULTS_DIR, 'windowHits/')        # 02: hit-level results, observed DIP only
SUMMARY_DIR <- paste0(RESULTS_DIR, 'summary/')           # 03
for (.d in c(COUNT_DIR, HIT_DIR, SUMMARY_DIR)) dir.create(.d, recursive = TRUE, showWarnings = FALSE)

# --- Prototype scope ---------------------------------------------------------------
TISSUE        <- .envPath('QTLSHUFFLE_TISSUE', 'Muscle_Skeletal')
N_GENES       <- .envNum('QTLSHUFFLE_N_GENES', 50)      # most rhythmic genes (smallest harmonic p)

# --- Analysis constants --------------------------------------------------------------
N_WINDOWS     <- 8          # non-overlapping 3 h DIP windows
MIN_WINDOW_N  <- 60         # windows with fewer donors are skipped (NA counts)
CIS_WINDOW    <- 1e6        # bp either side of the TSS
MAF_THRESHOLD <- 0.01       # minimum MAF across the tissue's donors
MIN_MAC       <- 10         # minimum minor allele count in each donor set tested
N_PEER        <- 10         # PEER factors, recomputed in every donor set (on all genes)
PEER_ITER     <- 1000       # PEER maximum iterations
NPERM_GENO    <- .envNum('QTLSHUFFLE_NPERM_GENO', 10000) # phenotype permutations per gene (tensorQTL default)
GENE_FDR      <- 0.05       # Storey q-value on beta-approximated p (GTEx)
QVALUE_LAMBDA <- 0.85       # q-value lambda used by GTEx
VANISH_P      <- 0.05       # 24 h reference lead-variant p above this ...
HET_P         <- 0.05       # ... and window vs reference heterogeneity p below this
N_PERM        <- .envNum('QTLSHUFFLE_N_PERM', 1000)    # DIP permutations (0 = observed)
PERM_BLOCK    <- .envNum('QTLSHUFFLE_PERM_BLOCK', 20)  # DIP permutations per SLURM task
SEED          <- 20261102

# PEER backend: 'peer' (the R package) is the analysis; 'pca_test_only' exists
# solely for the toy test where PEER is not installed.
PEER_BACKEND  <- .envPath('QTLSHUFFLE_PEER_BACKEND', 'peer')
# PEER is unmaintained and only builds under R 3.4 (bioconda r-peer 1.3), so it
# lives in its own conda environment, installed by the rhyh2 repo's
# Code/peer/install_peer.sh. Each fit runs there as a separate process.
PEER_RSCRIPT  <- .envPath('QTLSHUFFLE_PEER_RSCRIPT', '/home/uqschauq/software/peer_env/bin/Rscript')

.a <- commandArgs(FALSE); .f <- sub('^--file=', '', .a[grep('^--file=', .a)])
CODE_DIR <- if (length(.f)) dirname(normalizePath(.f)) else normalizePath('.')
source(file.path(CODE_DIR, '00_functions.R'))
