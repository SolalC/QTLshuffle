#!/bin/bash
# Runs the whole pipeline on simulated toy data (no real data) and checks the
# output. PEER is replaced by PCA here because the PEER package is not needed
# to test the plumbing; the analysis itself always uses PEER.
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
toy="$here/toy"
rm -rf "$toy"; mkdir -p "$toy"
cd "$here/.."

Rscript tests/make_toy_data.R "$toy"

export QTLSHUFFLE_ROOT="$toy/"
export QTLSHUFFLE_RESULTS="$toy/Results/"
export QTLSHUFFLE_GTEX_EXPR_DIR="$toy/expr/"
export QTLSHUFFLE_GTEX_COV_DIR="$toy/cov/"
export QTLSHUFFLE_GENOTYPE_PLINK="$toy/genotypes.rds"
export QTLSHUFFLE_DIP_FILE="$toy/DIP.RData"
export QTLSHUFFLE_HARMONIC_DIR="$toy/harmonic/"
export QTLSHUFFLE_N_TISSUES=2
export QTLSHUFFLE_N_PERM=39
export QTLSHUFFLE_PERM_BLOCK=20
export QTLSHUFFLE_PEER_BACKEND=pca_test_only

Rscript 01_prepareGenotypes.R
for i in 1 2; do Rscript 02_prepareTissues.R $i 2; done
for task in 1 2 3 4; do Rscript 03_windowScan.R $task 4; done
Rscript 04_summarise.R
Rscript tests/check_toy_results.R "$QTLSHUFFLE_RESULTS"
