#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=120G
#SBATCH --job-name=roll05_referenceRetest
#SBATCH --time=8:00:00
#SBATCH --partition=general
#SBATCH --constraint=epyc3
#SBATCH --account=a_imb_cpdg
#SBATCH --array 1-1000
#SBATCH -o ../logs/%x_%a.out
#SBATCH -e ../logs/%x_%a.err

# Submit after 04: sbatch --dependency=afterok:<04> 05_referenceRetest.sh
# Array 1-1000: task r = reference sample r, for all 24 centres (24 PEER fits in parallel; ~3 h, as 01 took 2 h 57 min for the same 24 fits).
# Time one task first: sbatch --array 1 05_referenceRetest.sh
cd "${SLURM_SUBMIT_DIR}"
module load r/4.4.2-heavy
Rscript 05_referenceRetest.R $SLURM_ARRAY_TASK_ID $SLURM_CPUS_PER_TASK
