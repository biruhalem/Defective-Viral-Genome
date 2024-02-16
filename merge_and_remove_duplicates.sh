#!/bin/bash

# Concatenate all files into one merged file
cat *_cb_100_r6.fa > merged_cb_100_r6.fa

# Extract unique sequences based on names and save to output file
awk 'BEGIN {RS=">"; FS="\n"} NR>1 {if (!seen[$1]++) printf ">%s", $0}' merged_cb_100_r6.fa > merge_unique_cb_100_r6.fa

