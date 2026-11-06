#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=120G
#SBATCH --job-name=roll04_windowCall
#SBATCH --time=12:00:00
#SBATCH --partition=general
#SBATCH --constraint=epyc3
#SBATCH --account=a_imb_cpdg

#SBATCH -o ../logs/%x.out
#SBATCH -e ../logs/%x.err

# Submit after every 03 task: sbatch --dependency=afterok:<03> 04_windowCall.sh
cd "${SLURM_SUBMIT_DIR}"
module load r/4.4.2-heavy
Rscript 04_windowCall.R $SLURM_CPUS_PER_TASK
