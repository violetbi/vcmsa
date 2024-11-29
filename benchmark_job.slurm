#!/bin/bash
#SBATCH --job-name=vcmsa_job   # Job name
#SBATCH --error=vcmsa_job.err  # Error log file
#SBATCH --ntasks=1             # Number of tasks (processes)
#SBATCH --cpus-per-task=1      # Number of CPU cores per task
#SBATCH --mem=128G              # Job memory request
#SBATCH --time=00:10:00        # Time limit hrs:min:sec
#SBATCH --mail-type=begin        # send email when job begins
#SBATCH --mail-type=end                # send email when job ends
#SBATCH --mail-user=xb4719@princeton.edu
#SBATCH --out=benchmark_job.out

# Load necessary modules (if needed)
module load python
# module load vcmsa

# Activate any environment (if needed)
source activate benchmark

# Function to compare alignments 
compare_alignments() {
    ref=$1
    alignments=("${@:2}")
    
    for aln in "${alignments[@]}"; do
        echo "Reference file: $ref"
        echo "Alignment file: $aln"

        for mode in tc; do
            echo "comparison $mode mode: "
            t_coffee -other_pg aln_compare -al1 $ref -al2 $aln -compare_mode $mode

            # capture output
            output=$(t_coffee -other_pg aln_compare -al1 $ref -al2 $aln -compare_mode $mode)
            
            all_score=$(echo "$output" | awk '/seq1/ {getline; print $4}')
            
            echo "[ALL] score for $aln: $all_score"
        done
    done
}

ref1="LDLa.ref"
alignments1=("LDLa.vie.20seqs.fasta.16layer.bert.aln" "LDLa.vie.20seqs.fasta.16layer.t5.aln" "MAFFT_GINSI_aligned_LDLa.vie.20seqs.fasta")

ref2="Asp_Glu_race_D.ref"
alignments2=("Asp_Glu_race_D.vie.20seqs.fasta.16layer.bert.aln" "Asp_Glu_race_D.vie.20seqs.fasta.16layer.t5.aln" "MAFFT_GINSI_aligned_Asp_Glu_race_D.vie.20seqs.fasta")

compare_alignments "$ref1" "${alignments1[@]}"
compare_alignments "$ref2" "${alignments2[@]}"

# # Set the reference file and alignment files
# ref=LDLa.ref
# alignments=("LDLa.vie.20seqs.fasta.16layer.bert.aln" "LDLa.vie.20seqs.fasta.16layer.t5.aln" "MAFFT_GINSI_aligned_LDLa.vie.20seqs.fasta")

# # Loop through each alignment and run comparisons
# for aln in "${alignments[@]}"; do
#     echo "Reference file: $ref"
#     echo "Alignment file: $aln"

#     for mode in tc; do
#     #for mode in sp tc column; do
#         echo "comparison $mode mode: "
#         t_coffee -other_pg aln_compare -al1 $ref -al2 $aln -compare_mode $mode

#         # capture output
#         output=$(t_coffee -other_pg aln_compare -al1 $ref -al2 $aln -compare_mode tc)
        
#         # Extract the [ALL] score using awk or grep
#         all_score=$(echo "$output" | awk '/seq1/ {getline; print $4}')
        
#         echo "[ALL] score for $aln: $all_score"

#     done
# done

# echo "conserverd poisitions between MSAs:"
# t_coffee -other_pg aln_compare -al1 $ref -al2 $aln -output_aln -output_aln_threshold 50
# t_coffee -other_pg aln_compare -al1 $ref -al2 $aln -output_aln -output_aln_threshold 50 -output_aln_modif x
