#!/bin/bash

model=../../prot_t5_xl_half_uniref50-enc
padding=0

st=0.7
max_iterations=100

# Loop over all .fasta files in the input directory
for input_fasta in ../../data/Fasta/*.fasta
do
    # Extract the base name of the fasta file (without the directory and extension)
    dataset_name=$(basename "$input_fasta" .fasta)
    output_dir=./output_$dataset_name
    mkdir -p $output_dir

    # layers='-16 -15 -14 -13 -12 -11 -10 -9 -8 -7 -6 -5 -4 -3 -2 -1'
    # output_file=$output_dir/layer-all.aln
    # echo "Running vcmsa for $dataset_name at all layers, output to $output_file"
    # vcmsa -m $model -l $layers -p $padding -st $st -mi $max_iterations -i $input_fasta -o $output_file --log ERROR

    # Loop over layers from -1 to -16
    for layer in {-1..-16}
    do
        output_file=$output_dir/layer${layer}.aln
        echo "Running vcmsa for $dataset_name at layer $layer, output to $output_file"
        vcmsa -m $model -l $layer -p $padding -st $st -mi $max_iterations -i $input_fasta -o $output_file --log ERROR
    done
done