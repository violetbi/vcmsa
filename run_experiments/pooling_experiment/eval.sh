#!/bin/bash

for input_fasta in ../../data/Fasta/*.fasta   #############################  Change this line
do

    # Extract the base name of the fasta file (without the directory and extension)
    dataset=$(basename "$input_fasta" .fasta)

    reference_file=../../data/Ref/${dataset}.ref

    # Output files for combined results
    sp_results_file=./results_${dataset}/all_sp_results.txt
    tc_results_file=./results_${dataset}/all_tc_results.txt

    mkdir -p ./results_${dataset}

    # Clear the result files if they already exist
    > $sp_results_file
    # > $column_results_file
    > $tc_results_file

    pooling=mean
    alignment_file=./output_${dataset}/$pooling.aln

    # Extract SP score
    sp_output=$(t_coffee -other_pg aln_compare -al1 $reference_file -al2 $alignment_file -compare_mode sp)
    sp_score=$(echo "$sp_output" | grep -v "seq1" | grep -v '*' | awk '{ print $4}' ORS="\t")
    echo "$pooling: $sp_score" >> $sp_results_file

    # Extract TC score
    tc_output=$(t_coffee -other_pg aln_compare -al1 $reference_file -al2 $alignment_file -compare_mode tc)
    tc_score=$(echo "$tc_output" | grep -v "seq1" | grep -v '*' | awk '{ print $4}' ORS="\t")
    echo "$pooling: $tc_score" >> $tc_results_file


    pooling=min
    alignment_file=./output_${dataset}/$pooling.aln

    # Extract SP score
    sp_output=$(t_coffee -other_pg aln_compare -al1 $reference_file -al2 $alignment_file -compare_mode sp)
    sp_score=$(echo "$sp_output" | grep -v "seq1" | grep -v '*' | awk '{ print $4}' ORS="\t")
    echo "$pooling: $sp_score" >> $sp_results_file

    # Extract TC score
    tc_output=$(t_coffee -other_pg aln_compare -al1 $reference_file -al2 $alignment_file -compare_mode tc)
    tc_score=$(echo "$tc_output" | grep -v "seq1" | grep -v '*' | awk '{ print $4}' ORS="\t")
    echo "$pooling: $tc_score" >> $tc_results_file


    pooling=max
    alignment_file=./output_${dataset}/$pooling.aln

    # Extract SP score
    sp_output=$(t_coffee -other_pg aln_compare -al1 $reference_file -al2 $alignment_file -compare_mode sp)
    sp_score=$(echo "$sp_output" | grep -v "seq1" | grep -v '*' | awk '{ print $4}' ORS="\t")
    echo "$pooling: $sp_score" >> $sp_results_file

    # Extract TC score
    tc_output=$(t_coffee -other_pg aln_compare -al1 $reference_file -al2 $alignment_file -compare_mode tc)
    tc_score=$(echo "$tc_output" | grep -v "seq1" | grep -v '*' | awk '{ print $4}' ORS="\t")
    echo "$pooling: $tc_score" >> $tc_results_file

done