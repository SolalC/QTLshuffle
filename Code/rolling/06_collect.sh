#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=60G
#SBATCH --job-name=roll06_collect
#SBATCH --time=2:00:00
#SBATCH --partition=general
#SBATCH --account=a_imb_cpdg
#SBATCH --constraint=epyc3
#SBATCH -o ../logs/%x.out
#SBATCH -e ../logs/%x.err

# Submit after every 05 task: sbatch --dependency=afterok:<05> 06_collect.sh
cd "${SLURM_SUBMIT_DIR}"
module load r/4.4.2-heavy
Rscript 06_collect.R
