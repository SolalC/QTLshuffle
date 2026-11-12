#!/bin/bash
# Runs the rolling-window pipeline on simulated toy data (no real data) and
# checks the output. Needs MatrixEQTL, BEDMatrix, genio and plink2 (PLINK2 in
# 00_config.R, or QTLSHUFFLE_PLINK2): run it on Bunya with
#   sbatch tests/run_toy_test_slurm.sh   (from Code/rolling)
# PEER is replaced by PCA here only to keep the test fast; the analysis uses PEER.
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
toy="$here/toy"
rm -rf "$toy"; mkdir -p "$toy"
cd "$here/.."
nCore=${NCORE:-4}

Rscript ../tests/make_toy_data.R "$toy"

export QTLSHUFFLE_ROOT="$toy/"
export QTLSHUFFLE_RESULTS="$toy/Results/"
export QTLSHUFFLE_GTEX_EXPR_DIR="$toy/expr/"
export QTLSHUFFLE_GTEX_COV_DIR="$toy/cov/"
export QTLSHUFFLE_GENOTYPE_PLINK="$toy/genotypes"   # toy plink file: 02a only splits it by chromosome
export QTLSHUFFLE_GENO_DIR="$toy/geno/"
export QTLSHUFFLE_DIP_FILE="$toy/DIP.RData"
export QTLSHUFFLE_TISSUE=ToyA
export QTLSHUFFLE_CHROMOSOMES=1,2
export QTLSHUFFLE_N_REF=40
export QTLSHUFFLE_N_NULL=10
export QTLSHUFFLE_PEER_BACKEND=pca_test_only

Rscript 01_prepareTissue.R $nCore
Rscript 02a_vcfToPlink.R 1 2000      # toy plink file exists: only the split by chromosome
for t in $(seq 1 48); do Rscript 03_windowScan.R $t 2>&1 | grep -E 'MatrixEQTL|skipped' || true; done
Rscript 04_windowCall.R $nCore
for r in $(seq 1 40); do Rscript 05_referenceRetest.R $r $nCore 2>&1 | grep -iE 'error' || true; done
Rscript 06_collect.R
for k in $(seq 1 10); do Rscript 07_nullScan.R $k $nCore 2>&1 | grep -iE 'error|centres scanned' || true; done
for c in $(seq 1 24); do Rscript 08_nullCalibrate.R $c $nCore 2>&1 | grep -iE 'error|Check:|calibrate|skipped' || true; done
Rscript tests/check_toy_results.R "$QTLSHUFFLE_RESULTS"
