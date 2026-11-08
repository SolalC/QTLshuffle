# =============================================================================
# Rolling-window eQTL (R only): one tissue, all expressed genes.
# Paths and constants for every script in this folder; sources 00_functions.R.
# Every input path can be overridden by an environment variable prefixed
# QTLSHUFFLE_ (used by tests/run_toy_test.sh).
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(MatrixEQTL)
})

.envDir <- function(var, default) {
  d <- Sys.getenv(var, default)
  if (!endsWith(d, '/')) d <- paste0(d, '/')
  d
}
.envPath <- function(var, default) Sys.getenv(var, default)
.envNum  <- function(var, default) as.numeric(Sys.getenv(var, default))

ROOT        <- .envDir('QTLSHUFFLE_ROOT', '/scratch/user/uqschauq/QTLshuffle/')
RESULTS_DIR <- .envDir('QTLSHUFFLE_RESULTS', paste0(ROOT, 'Results/rolling/'))

# --- Inputs (same files as the cd-eQTL analysis) ---------------------------------
CIRCADIAN_ROOT <- '/scratch/project_mnt/S0007/solal/circadian/'
# GTEx v10 normalised expression (inverse-normal transformed across all donors)
GTEX_EXPR_DIR  <- .envDir('QTLSHUFFLE_GTEX_EXPR_DIR', '/QRISdata/Q8106/Public/QTL/eQTL_expression_matrices/')
# GTEx v10 covariates; only sex, platform, pcr and PC1-5 are used (PEER is recomputed)
GTEX_COV_DIR   <- .envDir('QTLSHUFFLE_GTEX_COV_DIR', '/QRISdata/Q8106/Public/QTL/eQTL_covariates/')
# GTEx v9 WGS genotypes, all variants (953 donors; 634 GB)
GENOTYPE_VCF   <- .envPath('QTLSHUFFLE_GENOTYPE_VCF',
  '/QRISdata/Q8106/Controlled/Genotype/GTEx_Analysis_2021-02-11_v9_WholeGenomeSeq_953Indiv.vcf.gz')
# The only genotype file this pipeline writes (the plink file of 02a) goes
# here, not to RESULTS_DIR (controlled-access data). Later steps read it
# directly with BEDMatrix.
GENO_DIR       <- .envDir('QTLSHUFFLE_GENO_DIR', '/QRISdata/Q5776/QTLshuffle/')
# plink2 binary (no Bunya module; same as the rhyh2 pipeline)
PLINK2         <- .envPath('QTLSHUFFLE_PLINK2', '/home/uqschauq/software/plink2')
# Donor inferred phase (CHIRAL): R object `phi`, radians
DIP_FILE       <- .envPath('QTLSHUFFLE_DIP_FILE', paste0(CIRCADIAN_ROOT, 'Data/CHIRAL/DIP.RData'))
# Annotation for figures only (07): gene symbols (GENCODE v50 GTF, as in rhyh2)
# and rsIDs (GTEx v9 WGS lookup table, rs_id_dbSNP155_GRCh38p13, as in cd-eqtl)
GENE_GTF       <- .envPath('QTLSHUFFLE_GENE_GTF', '/QRISdata/Q9564/GENCODE/v50/GRCh38/gencode.v50.annotation.gtf.gz')
VARIANT_LOOKUP <- .envPath('QTLSHUFFLE_VARIANT_LOOKUP', paste0('/QRISdata/Q8106/Controlled/Genotype/support_files/',
                           'GTEx_Analysis_2021-02-11_v9_WholeGenomeSeq_953Indiv.lookup_table.txt.gz'))

# --- Scope -------------------------------------------------------------------------
TISSUE         <- .envPath('QTLSHUFFLE_TISSUE', 'Muscle_Skeletal')
# plink .bed/.bim/.fam made from GENOTYPE_VCF by 02a (tissue donors, PASS,
# biallelic, chr1-22, MAF >= MAF_THRESHOLD among the tissue donors), and one
# file per chromosome (02a) for plink2 --clump in 04
GENOTYPE_PLINK <- .envPath('QTLSHUFFLE_GENOTYPE_PLINK', paste0(GENO_DIR, 'plink/GTEx_v9_WGS_', TISSUE))
PLINK_CHR      <- \(chr) paste0(GENO_DIR, 'plink/', TISSUE, '_chr/', basename(GENOTYPE_PLINK), '_chr', chr)
CHROMOSOMES    <- as.integer(str_split_1(.envPath('QTLSHUFFLE_CHROMOSOMES', paste(1:22, collapse = ',')), ','))

# --- Windows -------------------------------------------------------------------------
WINDOW_WIDTH_H <- .envNum('QTLSHUFFLE_WINDOW_WIDTH_H', 3)  # window width (h of DIP)
CENTRES_H      <- 0:23      # window centres (h); windows wrap around midnight
MIN_WINDOW_N   <- .envNum('QTLSHUFFLE_MIN_WINDOW_N', 60)  # centres with fewer donors are skipped

# --- Association (MatrixEQTL, nominal threshold; no permutations) ---------------------
CIS_WINDOW     <- 1e6       # bp either side of the TSS
MAF_THRESHOLD  <- 0.01      # minimum MAF across the tissue's donors (02a)
MIN_MAC        <- 10        # minimum minor allele count in each donor set tested
PAIR_FDR       <- 0.05      # significant pair: Benjamini-Hochberg q < PAIR_FDR over all cis
                            # tests of the window (MatrixEQTL's FDR, pooled across chromosomes)
PAIR_P_STORE   <- 0.01      # 03 stores pairs with p below this (all pairs with q < PAIR_FDR
                            # have far smaller p; 04 checks)
CLUMP_R2       <- 0.01      # independent signals: plink2 --clump, r2 threshold
CLUMP_KB       <- 2000      # plink2 --clump window: the whole +/- 1 Mb cis window

# PEER factors as GTEx, by the number of donors in the set
peerK <- function(n) case_when(n < 150 ~ 15, n < 250 ~ 30, n < 350 ~ 45, TRUE ~ 60)
PEER_ITER      <- 1000      # PEER maximum iterations

# --- 24 h reference samples ------------------------------------------------------------
N_REF          <- .envNum('QTLSHUFFLE_N_REF', 1000)     # random 24 h samples per window centre (one SLURM task each)
SEED           <- 20261105

# --- Selection-bias control (07, 08) ------------------------------------------------
# References 1..N_NULL of each centre are treated as windows (scan, BH, clump:
# pseudo-windows); their signals are tested against the other N_NULL - 1, and
# the window signals against all N_NULL, with the same MatrixEQTL model.
N_NULL         <- .envNum('QTLSHUFFLE_N_NULL', 100)
# Bins of window |z| (|beta / se|) within which window signals are calibrated
# against pseudo-window signals (the bias depends on |z|)
Z_BINS         <- c(0, 5, 6, 7, 8, 10, 15, Inf)
TIME_FDR       <- 0.05      # time dependent: BH q (calibrated p, within window) < TIME_FDR

# PEER backend: 'peer' (the R package) is the analysis; 'pca_test_only' exists
# solely for the toy test where PEER is not installed.
PEER_BACKEND   <- .envPath('QTLSHUFFLE_PEER_BACKEND', 'peer')
# PEER only builds under R 3.4 (bioconda r-peer 1.3), so it lives in its own conda
# environment (rhyh2 repo, Code/peer/install_peer.sh); each fit is a separate process.
PEER_RSCRIPT   <- .envPath('QTLSHUFFLE_PEER_RSCRIPT', '/home/uqschauq/software/peer_env/bin/Rscript')

# --- Results ---------------------------------------------------------------------------
DATA_DIR  <- paste0(RESULTS_DIR, 'data/')                            # 01 (window-free)
WIDTH_DIR <- paste0(RESULTS_DIR, 'w', WINDOW_WIDTH_H, 'h/')          # everything per width
TISSUE_FILE  <- paste0(DATA_DIR, 'tissue.rds')                       # 01
WINDOW_FILE  <- paste0(WIDTH_DIR, 'windowCovariates.rds')            # 01
SCAN_DIR  <- paste0(WIDTH_DIR, 'scan/')                              # 03
CALL_DIR  <- paste0(WIDTH_DIR, 'calls/')                             # 04
REF_DIR   <- paste0(WIDTH_DIR, 'reference/')                         # 05: re-test results
REF_COV_DIR <- paste0(WIDTH_DIR, 'referenceCovariates/')            # 05: donors, PEER, covariates
OUT_DIR   <- paste0(WIDTH_DIR, 'results/')                           # 06
NULL_SCAN_DIR <- paste0(WIDTH_DIR, 'null/scan/')                     # 07: pseudo-window signals
NULL_CAL_DIR  <- paste0(WIDTH_DIR, 'null/calibration/')              # 08: re-tests, calibration
for (.d in c(DATA_DIR, SCAN_DIR, CALL_DIR, REF_DIR, REF_COV_DIR, OUT_DIR, NULL_SCAN_DIR, NULL_CAL_DIR)) dir.create(.d, recursive = TRUE, showWarnings = FALSE)

.a <- commandArgs(FALSE); .f <- sub('^--file=', '', .a[grep('^--file=', .a)])
CODE_DIR <- Sys.getenv('QTLSHUFFLE_CODE_DIR', if (length(.f)) dirname(normalizePath(.f)) else normalizePath('.'))
source(file.path(CODE_DIR, '00_functions.R'))
