#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=60G
#SBATCH --job-name=02_peerCovariates
#SBATCH --time=48:00:00
#SBATCH --partition=general
#SBATCH --array 1-510
#SBATCH --constraint=epyc3
#SBATCH --account=a_imb_cpdg
#SBATCH -o logs/%x_%a.out
#SBATCH -e logs/%x_%a.err

# Submit from this folder after 01 has finished: sbatch 02_peerCovariates.sh
# Array size = N_TISSUES x ceiling((N_PERM + 1) / PERM_BLOCK) = 10 x 51 = 510.
# Each task fits 2 x N_WINDOWS PEER models per permutation (320 for a block of
# 20). Time one task first (sbatch --array 1 02_peerCovariates.sh; the log
# prints minutes per permutation) and adjust --time and PERM_BLOCK.
cd "${SLURM_SUBMIT_DIR}"
module load r/4.4.2-heavy
Rscript 02_peerCovariates.R $SLURM_ARRAY_TASK_ID $SLURM_CPUS_PER_TASK
