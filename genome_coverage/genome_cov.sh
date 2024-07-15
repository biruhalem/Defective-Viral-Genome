#Coverage analysis for minus and plus strand reads for calculating the ratio of cbDVGs with viral genome 


#######################
#1.Filtering the minus and the plus tags from the .bam files
#######################

#for single end reads aligned .bam files 

#Ref https://samtools.github.io/hts-specs/SAMv1.pdf and  # https://www.biostars.org/p/196746/ #https://en.wikipedia.org/wiki/SAM_(file_format)

#!/bin/bash

# Iterate over each .bam file in the current directory
for bam_file in *.bam; do
    # Extract the file name without extension
    filename=$(basename -- "$bam_file")
    filename_no_ext="${filename%.bam}"

    # Filter plus strand reads
    samtools view -b -F 16 "$bam_file" > "${filename_no_ext}_plus.bam"

    # Filter minus strand reads
    samtools view -b -f 16 "$bam_file" > "${filename_no_ext}_minus.bam"
done

#################
#2. Calculate the coverage of each flag using samtools
#################
#!/bin/bash
# Iterate over .bam files
for bam_file in *_minus.bam; do
    # Define the output filename
    output_file="${bam_file%.bam}_coverage.txt"

for bam_file in *_plus.bam; do
    # Define the output filename
    output_file="${bam_file%.bam}_coverage.txt"

    # Calculate coverage depth using samtools depth
    samtools depth -a "$bam_file" > "$output_file"
done

######################################
#3.  Average coverage of the 1st 10K genomic region
######################################

import os

# Define a function to calculate the average of the first 10,000 values based on column 2 ID
def calculate_average(file_path):
    data = {}  # Dictionary to store data based on column 2 ID
    count = 0
    with open(file_path, 'r') as file:
        for line in file:
            parts = line.strip().split('\t')
            if len(parts) == 3:
                id_val = int(parts[1])
                if id_val <= 10000:
                    count += 1
                    if parts[0] in data:
                        data[parts[0]].append(int(parts[2]))
                    else:
                        data[parts[0]] = [int(parts[2])]
                else:
                    break
    averages = {}
    for key, value in data.items():
        averages[key] = round(sum(value) / len(value))  # Round the average value
    return averages, count

# Directory containing the *_Coverage.txt files
directory = os.getcwd()

# Output file to store the average values
output_file = 'plus_genomecoverage10k.txt'

# Open the output file in write mode
with open(output_file, 'w') as output:
    # Write headers for the table
    output.write("File\tSample\tAverage Value\n")

    # Iterate through each file in the directory
    for filename in os.listdir(directory):
        if filename.endswith('plus_coverage.txt'):
            file_path = os.path.join(directory, filename)
            averages, count = calculate_average(file_path)
            # Write the file name
            output.write(f"{filename}\n")
            # Write the average values for the file
            for key, value in averages.items():
                output.write(f"\t{key}\t{value}\n")
            output.write("\n")  # Add a newline between files

