#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=120G
#SBATCH --job-name=roll01_prepareTissue
#SBATCH --time=12:00:00
#SBATCH --partition=general
#SBATCH --constraint=epyc3
#SBATCH --account=a_imb_cpdg

#SBATCH -o logs/%x.out
#SBATCH -e logs/%x.err

# Submit from this folder: sbatch 01_prepareTissue.sh
# 24 PEER fits (one per window centre) run in parallel, one per core.
cd "${SLURM_SUBMIT_DIR}"
module load r/4.4.2-heavy
Rscript 01_prepareTissue.R $SLURM_CPUS_PER_TASK
