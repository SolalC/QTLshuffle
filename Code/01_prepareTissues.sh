#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=50G
#SBATCH --job-name=01_prepareTissues
#SBATCH --time=2:00:00
#SBATCH --partition=general
#SBATCH --array 1-10
#SBATCH --constraint=epyc3
#SBATCH --account=a_imb_cpdg
#SBATCH -o logs/%x_%a.out
#SBATCH -e logs/%x_%a.err

# Submit from this folder: sbatch 01_prepareTissues.sh
# (array size = N_TISSUES in 00_config.R)
cd "${SLURM_SUBMIT_DIR}"
module load r/4.4.2-heavy
Rscript 01_prepareTissues.R $SLURM_ARRAY_TASK_ID
