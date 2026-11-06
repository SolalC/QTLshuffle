#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --job-name=roll03_windowScan
#SBATCH --time=3:00:00
#SBATCH --partition=general
#SBATCH --constraint=epyc3
#SBATCH --account=a_imb_cpdg
#SBATCH --array 1-528
#SBATCH -o ../logs/%x_%a.out
#SBATCH -e ../logs/%x_%a.err

# Submit after 01 and 02a: sbatch --dependency=afterok:<01>:<02a> 03_windowScan.sh
# Array 1-528 = 24 centres x 22 chromosomes; task k: centre (k-1) %/% 22, chromosome (k-1) %% 22 + 1.
# MatrixEQTL, one core. Time one task first (chr1): sbatch --array 1 03_windowScan.sh
cd "${SLURM_SUBMIT_DIR}"
module load r/4.4.2-heavy
Rscript 03_windowScan.R $SLURM_ARRAY_TASK_ID
