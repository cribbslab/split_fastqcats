#Aim: check LR RNA seq fastq for correct primer orientation and trim reads/split them when concatenation has occurred
##This script defines QC based on the following parameters:
##-ONT adapter sequences are often missing from read due to removal during basecalling so this is not a parameter in this pipeline
##-Lack of polyA/polyT tail can sometimes occur and so will not qualify read for 'bin', but is used to orientate reads (reverse complementing when polyT present)
##-PCR primers should enclose the read
##-Homotrimer UMIs should be within the PCR primer and reads
##-Indels or single base substitutions are common in ONT seq so this is allowed for

##The tool will produce two outputs: 'bin' folder with reads chucked out, 'processed' output with trimmed and correct reads
##The tool will also give basic stats for use of metrics

#Setup
import argparse
from Bio import SeqIO #For handling fastq file parsing and writing 
from Bio.SeqRecord import SeqRecord #For creating new sequences
from Bio.Seq import Seq #For working with biological seqs
import gzip #Compressed file handling
import regex #Fuzzy matching
import csv #Writing stats to csv
from collections import defaultdict #Creating a dictionary for stats file
from Bio import Align
import numpy as np

def smith_waterman_search(sequence, primer, max_mismatches=5):
    """
    Smith-Waterman algorithm for local sequence alignment.
    Returns start, end, score, and number of errors.
    """
    # Initialise the local pairwise aligner
    aligner = Align.PairwiseAligner()
    aligner.mode = 'local'
    
    # Scoring:
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -1
    aligner.extend_gap_score = -0.5
    
    matches = []
    primer_length = len(primer)
    
    # Minimum score threshold of 75% match
    min_score_threshold = 0.75 * primer_length * aligner.match_score  

    # Loop that scans the sequence in windows matching the primer length
    for i in range(len(sequence) - primer_length + 1):
        sub_seq = sequence[i:i + primer_length]
        alignments = aligner.align(sub_seq, primer)
        
        if not alignments:
            continue  # Skip this iteration
        alignment = alignments[0]
        score = alignment.score
        
        
        # Filter based on score (counts mismatches and keeps acceptable mismatch cases)
        if score >= min_score_threshold:
            mismatches = sum(1 for a, b in zip(sub_seq, primer) if a != b)
            if mismatches <= max_mismatches:
                matches.append({'start': i, 'end': i + primer_length, 'score': score, 'mismatches': mismatches})
                
    # Sort by mismatch rate (ascending) and start position (ascending)
    return sorted(matches, key=lambda x: (x['mismatches'], x['start']))

def find_best_primer_pairs(sequence, forward_primer, reverse_primer):
    """
    Finds the closest forward and reverse primer pairs in a sequence using Smith-Waterman.
    """
    # Find forward and reverse primer matches using above function
    forward_matches = smith_waterman_search(sequence, forward_primer)
    reverse_matches = smith_waterman_search(sequence, reverse_primer)
    
    paired_sequences = []
    used_reverse_indices = set()

    for f_match in forward_matches:
        best_pair = None
        best_distance = float('inf')
        best_mismatches = float('inf')

        for r_match in reverse_matches:
            # Only considers reverse primers that follow a forward primer
            if r_match['start'] > f_match['end'] and r_match['start'] not in used_reverse_indices:
                distance = r_match['start'] - f_match['end']
                total_mismatches = f_match['mismatches'] + r_match['mismatches']
                
                # Choses primer pairs based on shortest distance and least mismatched
                if distance < best_distance or (distance == best_distance and total_mismatches < best_mismatches):
                    best_distance = distance
                    best_mismatches = total_mismatches
                    best_pair = r_match
        
        if best_pair:
            # Saves best primer pair and relevant info
            paired_sequences.append({
                'trimmed_seq': sequence[f_match['end']:best_pair['start']],
                'forward': (f_match['start'], f_match['end']),
                'reverse': (best_pair['start'], best_pair['end']),
                'total_mismatches': f_match['mismatches'] + best_pair['mismatches']
            })
            used_reverse_indices.add(best_pair['start'])

    return paired_sequences


def split_reads(
    input_file, forward_primer, reverse_primer, index_dictionary,
    processed_output, lowqual_output, bin_output, stats_output):
    
    total_sequences = 0
    sequences_with_one_match = 0
    total_forward_reverse_matches = 0
    sequences_with_one_index = 0
    sequences_with_two_indexes_far = 0
    total_processed = 0
    total_binned = 0
    total_lowqual = 0
    processed_records = []
    lowqual_records = []
    binned_records = []
    
    with gzip.open(input_file, 'rt') as handle:
        records = SeqIO.parse(handle, 'fastq')
        
        for record in records:
            total_sequences += 1
            seq = str(record.seq)
            qual = record.letter_annotations["phred_quality"]

            matches = find_best_primer_pairs(seq, forward_primer, reverse_primer)
            match_count = len(matches)
            total_forward_reverse_matches += match_count
        
            if match_count == 0 or match_count > 10:
                binned_records.append(record)
                total_binned += 1
                continue

            if match_count == 1:
                sequences_with_one_match += 1

            for match in matches:
                trimmed_seq = match['trimmed_seq']
                trimmed_qual = record.letter_annotations["phred_quality"][:len(trimmed_seq)]  # Get the quality scores corresponding to the trimmed sequence
            
                forward_index_found, reverse_index_found = False, False
                forward_index_number, reverse_index_number = None, None

                # Search for indexes in the trimmed sequence
                for index, (start_index, _) in index_dictionary.items():
                    
                    # Find matches for the forward index (start_index)
                    f_matches = smith_waterman_search(trimmed_seq, start_index)
                    if f_matches:
                        best_f_match = f_matches[0]
                        f_start = best_f_match['start']
                        f_end = best_f_match['end']
                        f_score = best_f_match['score']
                        if f_score >= 0.83 * len(start_index) * 2:
                            forward_index_found = True
                            forward_index_number = index
                            break  # Exit the loop once a forward index is found

                for index, (_, end_index) in index_dictionary.items():
                    
                    # Find matches for the reverse index (end_index)
                    r_matches = smith_waterman_search(trimmed_seq, end_index)
                    if r_matches:
                        best_r_match = r_matches[0]
                        r_start = best_r_match['start']
                        r_end = best_r_match['end']
                        r_score = best_r_match['score']
                        if r_score >= 0.83 * len(end_index) * 2:
                            reverse_index_found = True
                            reverse_index_number = index
                            break  # Exit the loop once a reverse index is found


                # Check positions and classify based on the new conditions
                if forward_index_found and reverse_index_found:
                    # Check if forward primer is within 5 bases of the start (position 0 to 5)
                    f_within_range = f_start is not None and f_start <= 5
                    # Check if reverse primer is within 5 bases of the end (last 5 positions)
                    r_within_range = r_end is not None and r_end >= len(trimmed_seq) - 5

                    if f_within_range and r_within_range:
                        # Both indexes are within range
                        processed_records.append(record)
                        total_processed += 1
                    else:
                        # Both indexes are found but not within range
                        lowqual_records.append(record)
                        total_lowqual += 1
                        sequences_with_two_indexes_far += 1
            
                elif forward_index_found or reverse_index_found:
                    # Only one index is found
                    lowqual_records.append(record)
                    total_lowqual += 1
                    sequences_with_one_index += 1
            
                else:
                    # No indexes found, bin the record
                    binned_records.append(record)
                    total_binned += 1
                   
            
    with gzip.open(processed_output, 'wt') as processed_handle:
        SeqIO.write(processed_records, processed_handle, 'fastq')

    with gzip.open(lowqual_output, 'wt') as lowqual_handle:
        SeqIO.write(lowqual_records, lowqual_handle, 'fastq')

    with gzip.open(bin_output, 'wt') as bin_handle:
        SeqIO.write(binned_records, bin_handle, 'fastq')

    # Write statistics to CSV
    with open(stats_output, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(["Metric", "Value"])
        writer.writerow(["Total sequences", total_sequences])
        writer.writerow(["Sequences with one primer match", sequences_with_one_match])
        writer.writerow(["Average primer matches per sequence", total_forward_reverse_matches / total_sequences])
        writer.writerow(["Total primer matches", total_forward_reverse_matches])
        writer.writerow(["Sequences with one index", sequences_with_one_index])
        writer.writerow(["Sequences with two indexes >5bp away", sequences_with_two_indexes_far])
        writer.writerow(["Processed sequences", total_processed])
        writer.writerow(["Low-quality sequences", total_lowqual])
        writer.writerow(["Binned sequences", total_binned])

# Run the function
split_reads(input_file, forward_primer, reverse_primer, indexes, processed_output, lowqual_output, bin_output, stats_output)
