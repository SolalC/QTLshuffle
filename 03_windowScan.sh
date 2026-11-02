#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=200G
#SBATCH --job-name=03_windowScan
#SBATCH --time=48:00:00
#SBATCH --partition=general
#SBATCH --array 1-510
#SBATCH --constraint=epyc3
#SBATCH --account=a_imb_cpdg
#SBATCH -o logs/%x_%a.out
#SBATCH -e logs/%x_%a.err

# Submit from this folder after 02 has finished: sbatch 03_windowScan.sh
# Array size = N_TISSUES x ceiling((N_PERM + 1) / PERM_BLOCK) = 10 x 51 = 510.
# Time one task first (sbatch --array 1 03_windowScan.sh; the log prints the
# minutes per permutation) and adjust --time, PERM_BLOCK and the array size.
# Memory: the tissue's genotypes (~20-30 GB) are shared by the forked workers.
cd "${SLURM_SUBMIT_DIR}"
module load r/4.4.0-combo-EPYC3-only
Rscript 03_windowScan.R $SLURM_ARRAY_TASK_ID $SLURM_CPUS_PER_TASK
