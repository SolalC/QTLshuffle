#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --job-name=roll02a_vcfToPlink
#SBATCH --time=24:00:00
#SBATCH --partition=general
#SBATCH --constraint=epyc3
#SBATCH --account=a_imb_cpdg
#SBATCH -o ../logs/%x.out
#SBATCH -e ../logs/%x.err

# Submit from this folder: sbatch 02a_vcfToPlink.sh   (no dependency on 01)
# Reads the whole 634 GB VCF once (skipped if the plink file exists), then
# splits it by chromosome. plink2: PLINK2 in 00_config.R.
cd "${SLURM_SUBMIT_DIR}"
module load r/4.4.2-heavy
Rscript 02a_vcfToPlink.R $SLURM_CPUS_PER_TASK 56000
