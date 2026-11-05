#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=60G
#SBATCH --job-name=roll02b_prepareGenotypes
#SBATCH --time=4:00:00
#SBATCH --partition=general
#SBATCH --constraint=epyc3
#SBATCH --account=a_imb_cpdg
#SBATCH --array 1-22
#SBATCH -o ../logs/%x_%a.out
#SBATCH -e ../logs/%x_%a.err

# Submit from this folder after 02a (not 01): sbatch --dependency=afterok:<02a> 02b_prepareGenotypes.sh
# One task per chromosome.
cd "${SLURM_SUBMIT_DIR}"
module load r/4.4.2-heavy
Rscript 02b_prepareGenotypes.R $SLURM_ARRAY_TASK_ID
