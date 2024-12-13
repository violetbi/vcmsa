model=prot_t5_xl_uniref50
layers='-16 -15 -14 -13 -12 -11 -10 -9 -8 -7 -6 -5 -4 -3 -2 -1'
padding=0

st=0.7
max_iterations=100

fasta=./example_files/LDLa.vie.20seqs.fasta

output_dir=./output

vcmsa  -m $model -l $layers -p $padding -st $st -mi $max_iterations -i $fasta -o $output_dir/alignment.aln --log INFO