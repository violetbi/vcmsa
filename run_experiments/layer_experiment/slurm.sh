#!/bin/bash
#SBATCH --job-name=vcmsa
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20G
#SBATCH --gres=gpu:1
#SBATCH --time=0:10:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=wz1411@princeton.edu

source ~/.bashrc
source /scratch/gpfs/wz1411/miniconda3/etc/profile.d/conda.sh

conda activate vcmsa_env

bash ./layer_experiment.sh