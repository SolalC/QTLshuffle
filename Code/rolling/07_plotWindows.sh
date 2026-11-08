#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=60G
#SBATCH --job-name=roll07_plot
#SBATCH --time=1:00:00
#SBATCH --partition=general
#SBATCH --account=a_imb_cpdg
#SBATCH --constraint=epyc3
#SBATCH -o ../logs/%x.out
#SBATCH -e ../logs/%x.err

# Run after 06: sbatch 07_plotWindows.sh [p=0.05] [test=pEmp|pEmpLarger] [independent=c,gene,var] [dependent=c,gene,var]
cd "${SLURM_SUBMIT_DIR}"
module load r/4.4.2-heavy
Rscript 07_plotWindows.R "$@"
