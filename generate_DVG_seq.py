#!/usr/bin/env python

from Bio import SeqIO

# Input files
sequence_file = "merge_unique_cb_100_r6_neg.fa"
lengths_file = "lengths_neg.txt"
output_file = "dvg_small_gen.fasta"

# Read the input sequence
records = list(SeqIO.parse(sequence_file, "fasta"))

# Open the output file
with open(output_file, 'w') as output_file:
    # Read the lengths from lengths_file
    with open(lengths_file, 'r') as lengths_file:
        for i, line in enumerate(lengths_file):
            # Get the corresponding record
            if i < len(records):
                seq_name = records[i].id
                length = int(line.strip())
                output_sequence = str(records[i].seq)[:length]
                output_file.write(f">{seq_name}_{length}\n{output_sequence}\n")

print("Done.")

