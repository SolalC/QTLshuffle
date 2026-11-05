#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --job-name=roll03_windowScan
#SBATCH --time=48:00:00
#SBATCH --partition=general
#SBATCH --constraint=epyc3
#SBATCH --account=a_imb_cpdg
#SBATCH --array 1-528
#SBATCH -o logs/%x_%a.out
#SBATCH -e logs/%x_%a.err

# Submit after 01 and 02: sbatch --dependency=afterok:<01>:<02> 03_windowScan.sh
# Array 1-528 = 24 centres x 22 chromosomes; task k: centre (k-1) %/% 22, chromosome (k-1) %% 22 + 1.
# Time one task first (chr1 is the largest): sbatch --array 1 03_windowScan.sh
cd "${SLURM_SUBMIT_DIR}"
module load r/4.4.2-heavy
Rscript 03_windowScan.R $SLURM_ARRAY_TASK_ID $SLURM_CPUS_PER_TASK
