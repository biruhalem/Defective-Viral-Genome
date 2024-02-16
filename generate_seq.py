#!/usr/bin/env python

from Bio import SeqIO

# Input files
sequence_file = "MVSSPE_G01.fasta"
lengths_file = "lengths.txt"
output_file = "MVSSPE_G01_small_seq.fasta"

# Read the lengths from lengths_file
with open(lengths_file, 'r') as lengths_file:
    lengths = [int(line.strip()) for line in lengths_file]

# Open the output file
with open(output_file, 'w') as output_file:
    # Read the input sequence
    for record in SeqIO.parse(sequence_file, "fasta"):
        seq_name = record.id
        sequence = str(record.seq)

        # Loop through the lengths
        for length in lengths:
            # Extract the first 'length' nucleotides
            output_sequence = sequence[:length]

            # Write the broken sequence to the output file
            output_file.write(f">{seq_name}_{length}\n{output_sequence}\n")

print("Done.")

