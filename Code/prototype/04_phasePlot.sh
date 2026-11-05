#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=50G
#SBATCH --job-name=proto04_phasePlot
#SBATCH --time=12:00:00
#SBATCH --partition=general
#SBATCH --constraint=epyc3
#SBATCH --account=a_imb_cpdg
#SBATCH -o ../logs/%x.out
#SBATCH -e ../logs/%x.err

# Submit from this folder after 01: sbatch 04_phasePlot.sh  (9 PEER fits)
cd "${SLURM_SUBMIT_DIR}"
module load r/4.4.2-heavy
Rscript 04_phasePlot.R
