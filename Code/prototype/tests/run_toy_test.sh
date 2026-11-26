#!/bin/bash
# Runs the prototype on simulated toy data (no real data) and checks the output.
# PEER is replaced by PCA here only because the PEER package is not installed
# where the test runs; the analysis itself uses PEER.
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
toy="$here/toy"
rm -rf "$toy"; mkdir -p "$toy"
cd "$here/.."

Rscript ../tests/make_toy_data.R "$toy"

export QTLSHUFFLE_ROOT="$toy/"
export QTLSHUFFLE_RESULTS="$toy/Results/"
export QTLSHUFFLE_GTEX_EXPR_DIR="$toy/expr/"
export QTLSHUFFLE_GTEX_COV_DIR="$toy/cov/"
export QTLSHUFFLE_GENOTYPE_PLINK="$toy/genotypes"
export QTLSHUFFLE_DIP_FILE="$toy/DIP.RData"
export QTLSHUFFLE_HARMONIC_DIR="$toy/harmonic/"
export QTLSHUFFLE_TISSUE=ToyA
export QTLSHUFFLE_N_PERM=39
export QTLSHUFFLE_PERM_BLOCK=20
export QTLSHUFFLE_NPERM_GENO=1000
export QTLSHUFFLE_PEER_BACKEND=pca_test_only

Rscript 01_prepareData.R
for b in 1 2; do Rscript 02_windowScan.R $b 4; done
Rscript 03_summarise.R
Rscript tests/check_toy_results.R "$QTLSHUFFLE_RESULTS"
