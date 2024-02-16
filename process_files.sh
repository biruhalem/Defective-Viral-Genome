#!/bin/bash

# List of input FASTA files
fasta_files=(*_fasta.fa)

# Loop through each FASTA file
for fasta_file in "${fasta_files[@]}"; do
    # Extract sequences based on positions in the corresponding text file
    base_name=$(basename "$fasta_file" .fa)
    output_file="${base_name}_cb_100_r6.fa"
    seqtk subseq "$fasta_file" "cb_dvg_r6_100.txt" > "$output_file"
done

