#!/bin/bash

# Usage: bash cat_files_by_index.sh --indexes "AAATTTGGGCCC TTTCCCAAAGGG ..."
# This script concatenates all files into merged ones.

INDEXES=$1

# Merge binned files
cat separate_samples.dir/*/*.binned_fastq.gz > merged_results.dir/all.binned_fastq.gz

# Merge lowqual files
cat separate_samples.dir/*/*.lowqual.fastq.gz > merged_results.dir/all.lowqual.fastq.gz

# Merge processed files per index
for index in $INDEXES; do
  cat separate_samples.dir/*/*_index_${index}.fastq.gz > merged_results.dir/all_index_${index}.fastq.gz
done

