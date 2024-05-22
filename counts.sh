#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 /path/to/parent_directory"
    exit 1
fi

# Set the parent directory from the argument
parent_dir="$1"

# Check if the provided directory exists
if [ ! -d "$parent_dir" ]; then
    echo "The directory $parent_dir does not exist."
    exit 1
fi

# Loop through each subdirectory in the parent directory
for sub_dir in "$parent_dir"/*/; do
    # Get the name of the subdirectory
    sub_dir_name=$(basename "$sub_dir")
    
    # Set the output file name based on the subdirectory name
    output_file="${parent_dir}/${sub_dir_name}.fastq"
    
    : > "$output_file" 
    
    # Merge all .fastq.gz files in the subdirectory into the output file
    for fastq_file in "$sub_dir"/*.fastq.gz; do
        zcat "$fastq_file" >> "$output_file"
    done
done

echo "Merging completed."

dir_name=$(basename "$parent_dir")

#counting function
perl ~/bin/count_ONT_barcodes.pl --seq_dir "$parent_dir" \
--barcodes ~/bin/barcode_lib.csv \
--ignore_missing_tag > "${parent_dir}/${dir_name}_counts.csv"
