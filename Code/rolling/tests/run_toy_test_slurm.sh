#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --job-name=rollToyTest
#SBATCH --time=2:00:00
#SBATCH --partition=general
#SBATCH --account=a_imb_cpdg
#SBATCH -o ../logs/%x.out
#SBATCH -e ../logs/%x.err

# Submit from Code/rolling: sbatch tests/run_toy_test.sh.sbatch
cd "${SLURM_SUBMIT_DIR}"
module load r/4.4.2-heavy
# Sourced, not run with bash: on Bunya, Rscript from the R module is only
# visible in this login shell, not in a child bash.
export NCORE=$SLURM_CPUS_PER_TASK
source tests/run_toy_test.sh
