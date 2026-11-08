#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=96G
#SBATCH --job-name=roll07_nullScan
#SBATCH --time=4:00:00
#SBATCH --partition=general
#SBATCH --constraint=epyc3
#SBATCH --account=a_imb_cpdg
#SBATCH --array 1-100
#SBATCH -o ../logs/%x_%a.out
#SBATCH -e ../logs/%x_%a.err

# Submit after 05 (needs referenceCovariates/ref<k>.rds): sbatch 07_nullScan.sh
# Array 1-100 = N_NULL: task k = reference sample k treated as a window at all 24 centres.
# Per task: 24 x 22 MatrixEQTL scans (centres in parallel) then plink2 clumping; no PEER.
# Time one task first: sbatch --array 1 07_nullScan.sh
cd "${SLURM_SUBMIT_DIR}"
module load r/4.4.2-heavy
Rscript 07_nullScan.R $SLURM_ARRAY_TASK_ID $SLURM_CPUS_PER_TASK
