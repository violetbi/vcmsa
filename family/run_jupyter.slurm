#!/bin/bash
#SBATCH --job-name=cath_mapping         # create a short name for your job
#SBATCH --nodes=1                      # node count
#SBATCH --ntasks=1                     # total number of tasks across all nodes
#SBATCH --cpus-per-task=2              # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=64G               # memory per cpu-core (adjust if needed)
#SBATCH --time=04:00:00                # total run time limit (adjust as needed)
#SBATCH --mail-type=begin        # send email when job begins
#SBATCH --mail-type=end                # send email when job ends
#SBATCH --mail-user=xb4719@princeton.edu
#SBATCH --out=run_jupyter.out


# # get tunneling info
# XDG_RUNTIME_DIR=""
# node=$(hostname -s)
# user=$(whoami)
# cluster="argo"
# port=8890

# print tunneling instructions jupyter-log
# echo -e "
# Command to create ssh tunnel:
# ssh -N -f -L ${port}:${node}:${port} ${user}@${cluster}.princeton.edu

# Use a Browser on your local machine to go to:
# localhost:${port}  (prefix w/ https:// if using password)
# "

# Run Jupyter
#jupyter-lab --no-browser --port=${port} --ip=${node}
# jupyter-notebook

# Load necessary modules and activate the conda environment
#module purge
#module load anaconda3/2024.6
# conda activate pytools-env

# source /Genomics/argo/users/xb4719/.conda/etc/profile.d/conda.sh  # Ensure this points to your Conda installation
# conda activate jupyter

# If you are running the Python script:
#python TRM_mapping.py

# Alternatively, if you want to run the notebook directly (no conversion to .py needed):
jupyter nbconvert --execute --inplace --to notebook cath_mapping.ipynb
#jupyter nbconvert --to script TRM_mapping.ipynb

#jupyter nbconvert --execute --inplace --to notebook check_adata.ipynb


