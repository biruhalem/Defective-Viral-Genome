#!/bin/bash

# Directories containing the fastq_pass files
directories=("SSPE_ONT_20240827T132846Z-001/Occ/fastq_pass"
             "SSPE_ONT_20240827T132846Z-001/Par/fastq_pass"
             "SSPE_ONT_20240827T132846Z-001/Temp/fastq_pass")

# Output directory for combined files
output_dir="combined_files300"
mkdir -p "$output_dir"

# Clean up any previous files in the output directory
rm -f "$output_dir"/*

# Step 1: Combine multiple runs of each sample into one
for dir in "${directories[@]}"; do
  uniquepathname=$(basename "$(dirname "$dir")")
  if ls "$dir"/*.fastq.gz 1> /dev/null 2>&1; then
    gunzip "$dir"/*.fastq.gz
  fi
  output_file="${output_dir}/merged_${uniquepathname}.fastq"
  rm -f "$output_file" # Remove old file if it exists
  cat "$dir"/*.fastq > "$output_file"
  echo "Combined files in $dir into $output_file"
done

# Step 2: Convert fastq reads from RNA to DNA
for fq in "${output_dir}"/*.fastq; do
  basename=$(basename "$fq" .fastq)
  output_file="${output_dir}/${basename}.dna.fastq"
  rm -f "$output_file" # Remove old file if it exists
  awk 'NR % 4 == 2 {gsub("U", "T"); print} NR % 4 != 2 {print}' "$fq" > "$output_file"
  echo "Converted $fq to $output_file"
done

# Step 3: Quality Trimming and Filtering
QUALITY_THRESHOLD=10
MIN_LENGTH=300
for fq in "${output_dir}"/*.dna.fastq; do
  basename=$(basename "$fq" .dna.fastq)
  output_file="${output_dir}/${basename}_quality.dna.fastq"
  rm -f "$output_file" # Remove old file if it exists
  NanoFilt -q "$QUALITY_THRESHOLD" -l "$MIN_LENGTH" < "$fq" > "$output_file"
  echo "Quality trimming and filtering complete for $fq, saved to $output_file"
done

# Step 4: Mapping to Reference Genome
GENOME_DIR="/mnt/c/Users/.../Desktop/SSPE_nanopore-20240718T172551Z-001"
GENOME_REF="MVSSPE_ancestor.fa.fasta"
for fq in "${output_dir}"/*_quality.dna.fastq; do
  basename=$(basename "$fq" _quality.dna.fastq)
  sam_output="${output_dir}/${basename}.sam"
  sorted_bam_output="${output_dir}/${basename}.sorted.bam"
  rm -f "$sam_output" "$sorted_bam_output" # Remove old files if they exist
  /mnt/c/Users/.../Desktop/SSPE_nanopore-20240718T172551Z-001/minimap2/minimap2 -ax map-ont "$GENOME_DIR/$GENOME_REF" "$fq" > "$sam_output"
  samtools view -bS -F 4 "$sam_output" | samtools sort -o "$sorted_bam_output"
  samtools index "$sorted_bam_output"
  rm "$sam_output"
  echo "Mapping and sorting complete for $fq, saved to $sorted_bam_output"
done

# Step 5: cbDVG Detection
COUNTS_FILE="${output_dir}/cbDVG_raw_result.txt"
rm -f "$COUNTS_FILE" # Remove old file if it exists
echo -e "cbDVG_ID_2\tcbDVG_ID_1\tbasename\treadname\tBreakpoint\tReinitiation\tstrand\truleofsix\tlength\tstem\tloop" > $COUNTS_FILE
for BAM_FILE in "${output_dir}"/*.sorted.bam; do
  BASENAME=$(basename "$BAM_FILE" .sorted.bam)
  FILTERED_READS="${output_dir}/${BASENAME}_filtered_reads.txt"
  rm -f "$FILTERED_READS" # Remove old file if it exists
  samtools view -F 4 "$BAM_FILE" | grep 'SA:Z' | awk '
  {
      primary_flag = $2
      sa_tag = $0
      split(sa_tag, arr, "SA:Z:")
      supplementary_info = arr[2]
      split(supplementary_info, sa_fields, ",")
      supplementary_flag = sa_fields[3]

      primary_strand = "+"
      if (and(primary_flag, 16)) {
          primary_strand = "-"
      }

      supplementary_strand = "+"
      if (supplementary_flag == "-") {
          supplementary_strand = "-"
      }

      if ((primary_strand == "+" && supplementary_strand == "-") || 
          (primary_strand == "-" && supplementary_strand == "+")) {
          print $0
      }
  }' > "$FILTERED_READS"

  if [[ ! -s "$FILTERED_READS" ]]; then
      continue
  fi

  while read -r READ; do
      READ_NAME=$(echo "$READ" | awk '{print $1}')
      REF_POS=$(echo "$READ" | awk '{print $4}')
      CIGAR=$(echo "$READ" | awk '{print $6}')
      SA_TAG=$(echo "$READ" | grep -o 'SA:Z:[^[:space:]]*')
      SUPPLEMENTARY_INFO=$(echo "$SA_TAG" | cut -d':' -f3)
      SUPPLEMENTARY_POS=$(echo "$SUPPLEMENTARY_INFO" | cut -d',' -f2)

      BREAKPOINT=$REF_POS
      REINITIATION=$SUPPLEMENTARY_POS
      CB_DVG_ID_1="cbDVG_${BREAKPOINT}_${REINITIATION}"

      STRAND=""
      if [[ $BREAKPOINT -gt $REINITIATION ]]; then
          STRAND="plus_strand"
          STEM=$((15894 - BREAKPOINT + 1))
          LOOP=$((BREAKPOINT - REINITIATION))
      else
          STRAND="minus_strand"
          STEM=$((15894 - REINITIATION + 1))
          LOOP=$((REINITIATION - BREAKPOINT))
      fi

      LENGTH=$((2 * STEM + LOOP))
      RULEOFSIX="No"
      if (( LENGTH % 6 == 0 )); then
          RULEOFSIX="Yes"
      fi

      CB_DVG_ID_2="cbDVG_${LENGTH}_${STEM}_${LOOP}"

      echo -e "$CB_DVG_ID_2\t$CB_DVG_ID_1\t$BASENAME\t$READ_NAME\t$BREAKPOINT\t$REINITIATION\t$STRAND\t$RULEOFSIX\t$LENGTH\t$STEM\t$LOOP" >> $COUNTS_FILE

  done < "$FILTERED_READS"

  rm "$FILTERED_READS"
done

# Step 6: Refining cbDVG Detection Results using Python
python3 <<EOF
import pandas as pd

# Load the input file
input_file = "${COUNTS_FILE}"
df = pd.read_csv(input_file, sep="\t")

# Drop the 'readname' column
df = df.drop(columns=["readname"])

# Create "Frequency_ID1" column (frequency of cbDVG_ID_1 in each basename)
df["Frequency_ID1"] = df.groupby(["cbDVG_ID_1", "basename"])["cbDVG_ID_1"].transform('count')

# Create "Frequency_ID2" column (frequency of cbDVG_ID_2 in each basename)
df["Frequency_ID2"] = df.groupby(["cbDVG_ID_2", "basename"])["cbDVG_ID_2"].transform('count')

# Assign True or False to the "both_strand" column
# True if both "minus_strand" and "plus_strand" are present for the same cbDVG_ID_2 in the same basename
df["both_strand"] = df.groupby(["cbDVG_ID_2", "basename"])["strand"].transform(lambda x: "minus_strand" in x.values and "plus_strand" in x.values)

# Ensure that the both_strand column contains only True or False values (not empty)
df["both_strand"] = df["both_strand"].fillna(False)

# Select rows where "ruleofsix" is "Yes"
df = df[df["ruleofsix"] == "Yes"]

# Select unique cbDVG_ID_1 in each basename
df_unique_ID1 = df.drop_duplicates(subset=["cbDVG_ID_1", "basename"])

# Select unique cbDVG_ID_2 in each basename based on strand selection criteria
def select_unique_cbDVG_ID_2(group):
    # If both_strand is True, select minus_strand
    if group["both_strand"].iloc[0]:
        minus_strand_group = group[group["strand"] == "minus_strand"]
        if not minus_strand_group.empty:
            return minus_strand_group.iloc[0]
    # If both_strand is False, select the cbDVG_ID_2 as is (no further filtering)
    return group.iloc[0]

# Apply the unique selection based on cbDVG_ID_2 and basename
df_unique_ID2 = df_unique_ID1.groupby(["cbDVG_ID_2", "basename"], as_index=False).apply(select_unique_cbDVG_ID_2)

# Sort by Frequency_ID2
df_unique_ID2 = df_unique_ID2.sort_values(by="Frequency_ID2", ascending=False)

# Save the final refined DataFrame to a file
output_file="${output_dir}/cbDVG_refined.txt"
df_unique_ID2.to_csv(output_file, sep="\t", index=False)
EOF

echo "Refinement process complete, saved to ${output_dir}/cbDVG_refined.txt"

