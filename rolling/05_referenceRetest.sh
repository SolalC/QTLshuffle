#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=60G
#SBATCH --job-name=roll05_referenceRetest
#SBATCH --time=48:00:00
#SBATCH --partition=general
#SBATCH --constraint=epyc3
#SBATCH --account=a_imb_cpdg
#SBATCH --array 1-480
#SBATCH -o logs/%x_%a.out
#SBATCH -e logs/%x_%a.err

# Submit after 04: sbatch --dependency=afterok:<04> 05_referenceRetest.sh
# Array 1-480 = 24 centres x 20 blocks of REF_BLOCK = 50 reference samples (N_REF = 1000).
# Each reference sample is one PEER fit; time one task first: sbatch --array 1 05_referenceRetest.sh
cd "${SLURM_SUBMIT_DIR}"
module load r/4.4.2-heavy
Rscript 05_referenceRetest.R $SLURM_ARRAY_TASK_ID $SLURM_CPUS_PER_TASK
