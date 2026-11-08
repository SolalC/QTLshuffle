#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=60G
#SBATCH --job-name=roll08_examples
#SBATCH --time=1:00:00
#SBATCH --partition=general
#SBATCH --account=a_imb_cpdg
#SBATCH --constraint=epyc3
#SBATCH -o ../logs/%x.out
#SBATCH -e ../logs/%x.err

# Run after 06: sbatch 08_plotExamples.sh [centre gene variant]
cd "${SLURM_SUBMIT_DIR}"
module load r/4.4.2-heavy
Rscript 08_plotExamples.R "$@"
