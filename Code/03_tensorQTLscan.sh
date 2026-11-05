#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --job-name=03_tensorQTLscan
#SBATCH --time=24:00:00
#SBATCH --partition=gpu_cuda
#SBATCH --gres=gpu:1
#SBATCH --array 1-510
#SBATCH --account=a_imb_cpdg
#SBATCH -o logs/%x_%a.out
#SBATCH -e logs/%x_%a.err

# Python environment with tensorQTL (and PyTorch built for CUDA):
TQ_ENV=/scratch/project_mnt/S0007/solal/envs/tensorqtl
# Configuration written by 01_prepareTissues.R (paths and constants from 00_config.R):
CONFIG=/scratch/user/uqschauq/QTLshuffle/Results/main/config.json

# Submit from this folder after the matching 02 tasks: sbatch 03_tensorQTLscan.sh
# Same task numbering as 02 (N_TISSUES x blocks = 510).
# TODO before the first run: check the Bunya GPU partition and --gres syntax,
# then time one task (sbatch --array 1 03_tensorQTLscan.sh; the log prints
# minutes per permutation). Each task runs 8 windows x PERM_BLOCK map_cis calls.
# Memory: the full genotype matrix is loaded once per task.
cd "${SLURM_SUBMIT_DIR}"
source "${TQ_ENV}/bin/activate"
python 03_tensorQTLscan.py "${CONFIG}" $SLURM_ARRAY_TASK_ID
