#!/bin/bash

# Usage: bash cat_files_by_index.sh --indexes "AAATTTGGGCCC TTTCCCAAAGGG ..."
# This script concatenates all files into merged ones.

INDEXES=$1
sample=$2

# Print values for debugging (optional)
echo "Using indexes: $INDEXES"
echo "Sample: $sample"


# Merge binned files
cat separate_samples.dir/${sample}*/*.binned_fastq.gz > merged_results.dir/${sample}.all.binned_fastq.gz

# Merge lowqual files
cat separate_samples.dir/${sample}*/*.lowqual.fastq.gz > merged_results.dir/${sample}.all.lowqual.fastq.gz

# Merge processed files per index
# Loop through each index in the INDEXES variable
IFS=' ' read -r -a index_array <<< "$INDEXES"  # Split the INDEXES string into an array

for index in "${index_array[@]}"; do
  # Check if any files match the pattern
  if ls separate_samples.dir/${sample}*/*_index_${index}.fastq.gz 1> /dev/null 2>&1; then
    # If files exist, concatenate them and output to the merged file
    cat separate_samples.dir/${sample}.*/*_index_${index}.fastq.gz > merged_results.dir/${sample}.all_index_${index}.fastq.gz
    echo "Files for index $index merged successfully."
  else
    # If no files exist for this index, print a message and skip
    echo "No files found for index $index. Skipping."
  fi
done


#for index in $INDEXES; do
#  cat separate_samples.dir/${sample}.*/*_index_${index}.fastq.gz > merged_results.dir/${sample}.all_index_${index}.fastq.gz
#done
