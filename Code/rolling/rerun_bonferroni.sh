#!/bin/bash
# Short rerun: re-call an existing run with the default rule (bonferroni_lead:
# Bonferroni p x nTests < 0.05, one lead pair per eGene), re-using its 03 scans
# and its 05 reference re-tests (no new scan, PEER fit or re-test):
# 04 -> 06 -> 09, chained with SLURM dependencies.
#
# The 05 output of a run made before the rule folders existed sits in
# w<W>h[_k<K>]/reference/ (BH + clumping); it is linked from
# bonferroni_lead/reference. This is valid because every bonferroni_lead signal
# was also a BH + clump signal (Bonferroni pairs pass BH; a gene's lead pair is
# its first clump index variant), so it was re-tested already.
#
# Usage (from Code/rolling):  bash rerun_bonferroni.sh <window width h> [fixed PEER k]
#   bash rerun_bonferroni.sh 3          # w3h/
#   bash rerun_bonferroni.sh 9 15       # w9h_k15/
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
W=${1:?window width in h}; K=${2:-}

# --- Paths to adjust ------------------------------------------------------------
RESULTS=${QTLSHUFFLE_RESULTS:-/scratch/user/uqschauq/QTLshuffle/Results/rolling/}
# ----------------------------------------------------------------------------------

export QTLSHUFFLE_WINDOW_WIDTH_H=$W
if [ -n "$K" ]; then export QTLSHUFFLE_PEER_K=$K; dir="${RESULTS%/}/w${W}h_k${K}"; else unset QTLSHUFFLE_PEER_K; dir="${RESULTS%/}/w${W}h"; fi
[ -d "$dir" ] || { echo "No results folder $dir"; exit 1; }

# 05 output of the earlier run (old layout, or already moved to BH_clump/)
ref=""
for d in "$dir/reference" "$dir/BH_clump/reference"; do [ -d "$d" ] && [ ! -L "$d" ] && { ref=$d; break; }; done
[ -n "$ref" ] || { echo "No 05 output in $dir/reference or $dir/BH_clump/reference"; exit 1; }
echo "$(ls "$ref" | grep -c '\.csv\.gz$' || true) reference files in $ref"
mkdir -p "$dir/bonferroni_lead"
[ -e "$dir/bonferroni_lead/reference" ] || ln -s "$ref" "$dir/bonferroni_lead/reference"

mkdir -p ../logs
for rule in bonferroni_lead; do
  export QTLSHUFFLE_CORRECTION=${rule%_*} QTLSHUFFLE_SIGNALS=${rule#*_}
  tag="w${W}h${K:+_k$K}_$rule"
  j4=$(sbatch --parsable --export=ALL -J "roll04_$tag" -o "../logs/roll04_$tag.out" -e "../logs/roll04_$tag.err" 04_windowCall.sh)
  j6=$(sbatch --parsable --export=ALL --dependency=afterok:$j4 -J "roll06_$tag" -o "../logs/roll06_$tag.out" -e "../logs/roll06_$tag.err" 06_collect.sh)
  j9=$(sbatch --parsable --export=ALL --dependency=afterok:$j6 -J "roll09_$tag" -o "../logs/roll09_$tag.out" -e "../logs/roll09_$tag.err" 09_plotWindows.sh)
  echo "$tag: 04 $j4 -> 06 $j6 -> 09 $j9"
done
