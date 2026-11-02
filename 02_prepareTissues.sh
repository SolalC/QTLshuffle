#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=100G
#SBATCH --job-name=02_prepareTissues
#SBATCH --time=12:00:00
#SBATCH --partition=general
#SBATCH --array 1-10
#SBATCH --constraint=epyc3
#SBATCH --account=a_imb_cpdg
#SBATCH -o logs/%x_%a.out
#SBATCH -e logs/%x_%a.err

# Submit from this folder after 01 has finished: sbatch 02_prepareTissues.sh
# (array size = N_TISSUES in 00_config.R)
cd "${SLURM_SUBMIT_DIR}"
module load r/4.4.0-combo-EPYC3-only
Rscript 02_prepareTissues.R $SLURM_ARRAY_TASK_ID $SLURM_CPUS_PER_TASK
