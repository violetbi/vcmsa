#!/bin/bash
#SBATCH --job-name=vcmsa_job   # Job name
#SBATCH --error=vcmsa_job.err  # Error log file
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1             # Number of tasks (processes)
#SBATCH --cpus-per-task=2      # Number of CPU cores per task
#SBATCH --mem=8G              # Job memory request
#SBATCH --gres=gpu:1             # number of gpus per node
#SBATCH --time=24:00:00        # Time limit hrs:min:sec
#SBATCH --mail-type=begin        # send email when job begins
#SBATCH --mail-type=end                # send email when job ends
#SBATCH --mail-user=xb4719@princeton.edu
#SBATCH --out=vcmsa_job.out

# Load necessary modules (if needed)
module load python
# module load vcmsa

# Activate any environment (if needed)
source activate vcmsa_env

# Set the parameters
#fasta=LDLa.vie.20seqs.fasta
#fasta=Asp_Glu_race_D.vie.20seqs.fasta
#fasta=QuanTest2/Test/Asp_Glu_race_D.vie.20seqs.fasta
layers='-16 -15 -14 -13 -12 -11 -10 -9 -8 -7 -6 -5 -4 -3 -2 -1'
model=prot_t5_xl_uniref50
suffix=16layer.t5.aln

#model=prot_bert_bfd
#suffix=.16layer.bert.aln

# Run the command
# vcmsa -i $fasta -o $fasta.$suffix --layers $layers -m $model --exclude --pca_plot
#vcmsa -i $fasta -o output/$(basename $fasta)$suffix --layers $layers -m $model --exclude --pca_plot

for fasta in QuanTest2/Test/*.fasta; do
    # Run the vcmsa command for each file
    vcmsa -i $fasta -o $fasta.$suffix --layers $layers -m $model --exclude --pca_plot
done