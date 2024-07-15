############################################
#Organizing the count.txt output data from DI-tector
###############################################
#!/bin/bash

# Specify the pattern to match the files
file_pattern="*counts.txt"

# Get a list of files matching the pattern
file_list=($file_pattern)

# Output merged file
merged_file="merged_output.txt"

# Loop through each file and apply the script
for filename in "${file_list[@]}"; do
    output_filename="${filename%.txt}_filtered.txt"  # Output filename based on the original filename

    cat "$filename" | grep -vE '^[[:space:]]*$|=' | awk 'NR==1 || (!/^=/ && !/^DVG'\''s/)' | awk 'BEGIN {OFS=FS="\t"} {print $1, $2, $3"_"$4, $5, $6, $7, $8, $9, $10}' | sed -e 's/ //g' > "$output_filename"

    echo "Processed $filename. Filtered content saved to $output_filename"
done

# Specify the pattern to match the filtered files
filtered_pattern="*filtered.txt"

# Output merged file
echo -e "DVG'stype\tLength\tBP_Pos_RI_Pos\tDelta_Positions\tRef\tCounts\t%_to_Virus" > "$merged_file"

# Loop through each filtered file and append based on BP_Pos_RI_Pos
for filtered_file in $filtered_pattern; do
    # Extract filename without extension
    filename=$(basename "$filtered_file" .txt)

    # Remove spaces between characters in the DVG's column
    sed -e 's/ //g' "$filtered_file" > "temp_filtered.txt"

    # Append lines based on BP_Pos_RI_Pos
    awk -v filename="$filename" 'NR>1 {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"filename}' "temp_filtered.txt" >> "$merged_file"

    # Remove the temporary file
    rm "temp_filtered.txt"
done

echo "Merged content saved to $merged_file"


