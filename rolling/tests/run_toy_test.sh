#!/bin/bash
# Runs the rolling-window pipeline on simulated toy data (no real data) and
# checks the output. PEER is replaced by PCA here only because the PEER package
# is not installed where the test runs; the analysis itself uses PEER.
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
toy="$here/toy"
rm -rf "$toy"; mkdir -p "$toy"
cd "$here/.."
nCore=${NCORE:-4}

Rscript ../tests/make_toy_data.R "$toy"

export QTLSHUFFLE_ROOT="$toy/"
export QTLSHUFFLE_RESULTS="$toy/Results/"
export QTLSHUFFLE_GTEX_EXPR_DIR="$toy/expr/"
export QTLSHUFFLE_GTEX_COV_DIR="$toy/cov/"
export QTLSHUFFLE_GENOTYPE_PLINK="$toy/genotypes"
export QTLSHUFFLE_DIP_FILE="$toy/DIP.RData"
export QTLSHUFFLE_TISSUE=ToyA
export QTLSHUFFLE_CHROMOSOMES=1,2
export QTLSHUFFLE_NPERM_GENO=1000
export QTLSHUFFLE_N_REF=40
export QTLSHUFFLE_REF_BLOCK=20
export QTLSHUFFLE_PEER_BACKEND=pca_test_only

Rscript 01_prepareTissue.R $nCore
for c in 1 2; do Rscript 02_prepareGenotypes.R $c; done
for t in $(seq 1 48); do Rscript 03_windowScan.R $t $nCore 2>&1 | grep -E 'tested|skipped' || true; done
Rscript 04_windowCall.R $nCore
for t in $(seq 1 48); do Rscript 05_referenceRetest.R $t $nCore 2>&1 | grep -iE 'error|no independent' || true; done
Rscript 06_collect.R
Rscript tests/check_toy_results.R "$QTLSHUFFLE_RESULTS"
