#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=120G
#SBATCH --job-name=roll08_nullCalibrate
#SBATCH --time=4:00:00
#SBATCH --partition=general
#SBATCH --constraint=epyc3
#SBATCH --account=a_imb_cpdg
#SBATCH --array 1-24
#SBATCH -o ../logs/%x_%a.out
#SBATCH -e ../logs/%x_%a.err

# Submit after every 07 task: sbatch --dependency=afterok:<07> 08_nullCalibrate.sh
# Array 1-24: task = window centre index. N_NULL references re-tested in parallel (MatrixEQTL on the
# window and pseudo-window pairs only), then calibration. Time one task first: sbatch --array 1 08_nullCalibrate.sh
cd "${SLURM_SUBMIT_DIR}"
module load r/4.4.2-heavy
Rscript 08_nullCalibrate.R $SLURM_ARRAY_TASK_ID $SLURM_CPUS_PER_TASK
