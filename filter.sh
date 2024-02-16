#!/bin/bash

# Specify the pattern to match the files
file_pattern="*counts.txt"

# Get a list of files matching the pattern
file_list=($file_pattern)

# Loop through each file and apply the script
for filename in "${file_list[@]}"; do
    output_filename="${filename%.txt}_filtered.txt"  # Output filename based on the original filename

    cat "$filename" | grep -vE '^[[:space:]]*$|=' | awk 'NR==1 || (!/^=/ && !/^DVG'\''s/)'|awk 'BEGIN {OFS=FS="\t"} {print $1, $2, $3"_"$4, $5, $6, $7, $8, $9, $10}'| sed -e 's/ //g' > "$output_filename"

    echo "Processed $filename. Filtered content saved to $output_filename"
done




