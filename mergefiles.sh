#!/bin/bash

# Specify the pattern to match the filtered files
filtered_pattern="*filtered.txt"

# Output merged file
merged_file="merged_output.txt"

# Create header with filenames
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

