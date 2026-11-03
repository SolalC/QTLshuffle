#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=60G
#SBATCH --job-name=proto02_windowScan
#SBATCH --time=24:00:00
#SBATCH --partition=general
#SBATCH --array 1-51
#SBATCH --constraint=epyc3
#SBATCH --account=a_imb_cpdg
#SBATCH -o logs/%x_%a.out
#SBATCH -e logs/%x_%a.err

# Submit from this folder after 01: sbatch 02_windowScan.sh
# Array size = ceiling((N_PERM + 1) / PERM_BLOCK) = 51 blocks of 20 DIP permutations.
# Each permutation fits 2 x 8 PEER models and 8 x N_GENES gene scans with
# NPERM_GENO permutations each. Time one task first (sbatch --array 1
# 02_windowScan.sh; the log prints minutes per permutation).
cd "${SLURM_SUBMIT_DIR}"
module load r/4.4.0-combo-EPYC3-only
Rscript 02_windowScan.R $SLURM_ARRAY_TASK_ID $SLURM_CPUS_PER_TASK
