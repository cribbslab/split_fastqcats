##Edits to Adam and Beth's main 'fastq_splitter_editted.py' script 
# - Reverse complement polyT sequences 
# - Remove UMI based QC processing - UMI can be dealt with by the tallytrin count pipeline
# - Filters reads without polyA
# . Error tolerance - -e parameter set number of mismatches/error tolerance - 0.3-0.4 seems to be a reasonable sensitivty for current test reads
# - Parallelised for efficient processing with chunk_size and num_workers arguments
# . Logging - logs to track reads and segments added
# . stats output - modified stats ouptput so reads and segments processed are recorded independent. 
# . Progress bar - added progress bar
# . Changes to find_best_primers function to remove the inner loop, added custom error tolerance. changed segment processing so all possible segments get classified as Valid/Invalid to be logged, and valid segments are more strictly defined. 

import os
import logging
import time
from tqdm import tqdm
import regex
import math
import sys
import argparse
import gzip
import parasail
from Bio import SeqIO, Align
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import csv
from typing import Dict, List, Tuple, Optional
from multiprocessing import Pool, cpu_count

class FastqSplitter:
    def __init__(self, forward_primer: str, reverse_primer: str, error: float, sw: int):
        """
        Initialize the FastQ splitter with primers and index dictionary.
        
        Args:
            forward_primer: Forward primer sequence
            reverse_primer: Reverse primer sequence
            index_dict: Dictionary mapping index numbers
            error: Mismatch allowance (float value)
        """
        
        self.forward_primer = forward_primer
        self.reverse_primer = reverse_primer
        self.error = error
        self.sw_version = sw
        self.aligner = Align.PairwiseAligner()
        self._setup_aligner()

    def _setup_aligner(self):
        """Configure the Smith-Waterman aligner with scoring parameters."""
        #Biopython PairwiseAligner arguments
        self.aligner.mode = 'local'
        self.aligner.match_score = 2
        self.aligner.mismatch_score = -1
        self.aligner.open_gap_score = -1
        self.aligner.extend_gap_score = -0.5
        #Parasail arguments
        self.match_score = 2
        self.mismatch_score = -1
        self.open_gap_score = 1
        self.extend_gap_score = 1
    
    #original single step sliding window Biopython aligner - very slow 10 reads/sec
    def smith_waterman_search(self, sequence: str, read_name: str, primer: str, error: float = 0.35) -> List[dict]:
        """
        Perform Smith-Waterman local alignment to find primer matches.
        
        Args:
            sequence: Input DNA sequence
            primer: Primer sequence to search for
            error: Mismatch allowance (float value)
            
        Returns:
            List of dictionaries containing match information
        """
        
        log_message(f"Performing Smith-Waterman search for Read={read_name} with primer={primer}", logging.DEBUG)
    
        
        matches = []
        primer_length = len(primer)
        
        # Determine min_score_threshold
        error = self.error
        if error < 1:
            min_score_threshold = (1 - error) * primer_length * self.aligner.match_score
            max_mismatches = math.ceil(error * primer_length)
        else:
            min_score_threshold = (1 - (error / primer_length)) * primer_length * self.aligner.match_score
            max_mismatches = int(error)  # Use error directly as max_mismatches
        
        ## uncomment below for testing puproses
        #min_score_threshold = 0.65 * primer_length * self.aligner.match_score  
        

        # Also search for reverse complement
        rev_comp_primer = str(Seq(primer).reverse_complement())
        
        for search_primer in [primer, rev_comp_primer]:
            
            for i in range(len(sequence) - primer_length + 1):
                sub_seq = sequence[i:min(i + primer_length, len(sequence))]
                alignments = self.aligner.align(sub_seq, search_primer)
                
                if not alignments:
                    continue
                    
                alignment = alignments[0]
                score = alignment.score
                
                if score >= min_score_threshold:
                    #max_mismatches = 10 #for testing purposes
                    mismatches = sum(1 for a, b in zip(sub_seq, search_primer) if a != b)
               
                    if mismatches <= max_mismatches:
                      
                        matches.append({
                            'start': i,
                            'end': i + primer_length,
                            'primer': search_primer,
                            'score': score,
                            'mismatches': mismatches,
                            'sequence': read_name,
                            'is_reverse': search_primer == rev_comp_primer
                        })
                    
        return sorted(matches, key=lambda x: (x['start'], -x['score'], x['mismatches']))

    ##Biopython aligner with 100bp windows and 50bp overlapping step - fast 125 reads per sec
    def smith_waterman_search2(self, sequence: str, read_name: str, primer: str, error: float = 0.35) -> List[dict]:
        """
        Perform Smith-Waterman local alignment to find primer matches.
        
        Args:
            sequence: Input DNA sequence
            primer: Primer sequence to search for
            error: Mismatch allowance (float value)
            
        Returns:
            List of dictionaries containing match information
        """
        
        log_message(f"Performing Smith-Waterman search for Read={read_name} with primer={primer}", logging.DEBUG)
    
        
        matches = []
        primer_length = len(primer)
        
        # Determine min_score_threshold
        error = self.error
        if error < 1:
            min_score_threshold = (1 - error) * primer_length * self.aligner.match_score
            max_mismatches = math.ceil(error * primer_length)
        else:
            min_score_threshold = (1 - (error / primer_length)) * primer_length * self.aligner.match_score
            max_mismatches = int(error)  # Use error directly as max_mismatches
        
        ## uncomment below for testing puproses
        #min_score_threshold = 0.65 * primer_length * self.aligner.match_score  
        

        # Also search for reverse complement
        rev_comp_primer = str(Seq(primer).reverse_complement())
        
        for search_primer in [primer, rev_comp_primer]:
            
            
            # set winow size for sliding search
            window_size = 100
            step_size = 50
            
            for i in range(0, len(sequence) - primer_length + 1, step_size):
                sub_seq = sequence[i:min(i + window_size, len(sequence))]
                alignments = self.aligner.align(sub_seq, search_primer)
                
                if not alignments:
                    continue
                  
                filtered_alignments = [alignment for alignment in alignments if alignment.score >= min_score_threshold]   
                
                if not filtered_alignments:
                    continue
                    
                alignment = filtered_alignments[0]
                score = alignment.score
                
                start_pos = alignment.coordinates[0][0] + i
                end_pos = min(len(sequence), start_pos+primer_length)
                target = sequence[start_pos:end_pos]
                query = search_primer[0:len(target)]
                
                if score >= min_score_threshold:
                    #max_mismatches = 10 #for testing purposes
                    mismatches = sum(1 for a, b in zip(target, query) if a != b)
               
                    if mismatches <= max_mismatches:
                      
                        matches.append({
                            'start': start_pos,
                            'end': end_pos,
                            'primer': search_primer,
                            'score': score,
                            'mismatches': mismatches,
                            'sequence': read_name,
                            'is_reverse': search_primer == rev_comp_primer
                        })
                    
        return sorted(matches, key=lambda x: (x['start'], -x['score'], x['mismatches']))

    #Biopython aligner without windows - simply take all alignments across full sequence - fast 250 reads per sec
    def smith_waterman_search3(self, sequence: str, read_name: str, primer: str, error: float = 0.35) -> List[dict]:
        """
        Perform Smith-Waterman local alignment to find primer matches.
        
        Args:
            sequence: Input DNA sequence
            primer: Primer sequence to search for
            error: Mismatch allowance (float value)
            
        Returns:
            List of dictionaries containing match information
        """
        
        log_message(f"Performing Smith-Waterman search for Read={read_name} with primer={primer}", logging.DEBUG)
    
        
        matches = []
        primer_length = len(primer)
        
        # Determine min_score_threshold
        error = self.error
        if error < 1:
            min_score_threshold = (1 - error) * primer_length * self.aligner.match_score
            max_mismatches = math.ceil(error * primer_length)
        else:
            min_score_threshold = (1 - (error / primer_length)) * primer_length * self.aligner.match_score
            max_mismatches = int(error)  # Use error directly as max_mismatches
        
        ## uncomment below for testing puproses
        #min_score_threshold = 0.65 * primer_length * self.aligner.match_score  
        

        # Also search for reverse complement
        rev_comp_primer = str(Seq(primer).reverse_complement())
        
        for search_primer in [primer, rev_comp_primer]:
        
            alignments = self.aligner.align(sequence, search_primer)
            
            if not alignments:
                #print("No alignment")
                continue
            
            filtered_alignments = [alignment for alignment in alignments if alignment.score >= min_score_threshold]   
            
            if not filtered_alignments:
                #print("No filtered alignment")
                continue
              
            for alignment in filtered_alignments:
                # Print the search primer (assume `search_primer` is defined elsewhere in the code)
                #print(search_primer)
                #print(alignment.target)
                start_pos = alignment.coordinates[0][0]
                end_pos = min(len(sequence), start_pos+primer_length)
                target = sequence[start_pos:end_pos]
                query = search_primer[0:len(target)]
                mismatches = sum(1 for a, b in zip(target, query) if a != b)
                if mismatches <= max_mismatches:
                  
                    matches.append({
                        'start': start_pos,
                        'end': end_pos,
                        'primer': search_primer,
                        'score': alignment.score,
                        'mismatches': mismatches,
                        'sequence': read_name,
                        'is_reverse': search_primer == rev_comp_primer
                    })
             
              
        return sorted(matches, key=lambda x: (x['start'], -x['score'], x['mismatches']))
    
    #Parasail aligner tool with 100bp window and 50bp overlapping step - very fast -- 1000 reads/sec
    def smith_waterman_search4(self, sequence: str, read_name: str, primer: str, error: float = 0.35) -> List[dict]:
        """
        Perform Smith-Waterman local alignment to find primer matches.
        
        Args:
            sequence: Input DNA sequence
            primer: Primer sequence to search for
            error: Mismatch allowance (float value)
            
        Returns:
            List of dictionaries containing match information
        """
        
        log_message(f"Performing Smith-Waterman search for Read={read_name} with primer={primer}", logging.DEBUG)
    
        
        matches = []
        primer_length = len(primer)
        seen_matches = set()
        
        # Determine min_score_threshold
        error = self.error
        if error < 1:
            min_score_threshold = (1 - error) * primer_length * self.aligner.match_score
            max_mismatches = math.floor(error * primer_length)
        else:
            min_score_threshold = (1 - (error / primer_length)) * primer_length * self.aligner.match_score
            max_mismatches = int(error)  # Use error directly as max_mismatches
        
        #min_score_threshold = 0.80 * primer_length * self.aligner.match_score  

        # Also search for reverse complement
        rev_comp_primer = str(Seq(primer).reverse_complement())
        scoring_matrix = parasail.matrix_create("ACGT", self.match_score, self.mismatch_score)
        
        
        # set winow size for sliding search
        window_size = 100
        step_size = 50
        
        for start in range(0, len(sequence) - primer_length + 1, step_size):
            window = sequence[start:min(start + window_size, len(sequence))]  # Extract window
        
            for search_primer in [primer, rev_comp_primer]:
            
                pattern_length = len(search_primer)
                
                # Perform Parasail alignment
                alignment = parasail.sw_trace_striped_16(
                    window, search_primer,
                    self.open_gap_score, # Gap opening penalty
                    self.extend_gap_score, # Gap extension penalty
                    scoring_matrix
                    #parasail.dnafull  # DNA-specific scoring matrix
                )
                
                # **Check if alignment score is valid**
                if alignment is None or alignment.score is None or alignment.traceback is None:
                    continue  # Skip this pattern if alignment failed
                
                # Ensure start position is valid
                start_position = max(0, alignment.end_query - pattern_length + 1 + start)
                end_position = min(len(sequence), start_position + pattern_length)
            
                # Create a unique identifier for each match
                match_key = (search_primer, start_position)
    
                # Avoid adding duplicates
                if match_key in seen_matches:
                    continue  # Skip duplicates
        
                # Otherwise, add to results
                seen_matches.add(match_key)
    
                if alignment.score >= min_score_threshold:
                    # Extract actual alignment details
                    aligned_query = alignment.traceback.query
                    aligned_pattern = alignment.traceback.ref
        
                    # Count mismatches and gaps
                    #mismatches = sum(1 for q, p in zip(aligned_query, aligned_pattern) if q != p and q != '-' and p != '-')
                    #gaps = sum(1 for q, p in zip(aligned_query, aligned_pattern) if q == '-' or p == '-')
                    mismatches = sum(1 for a, b in zip(aligned_query, aligned_pattern) if a != b)
                    #gaps = sum(1 for q, p in zip(aligned_query, aligned_pattern) if q == '-' or p == '-')
                    #mismatches = self.count_mismatches(aligned_query, aligned_pattern)
                    gaps = aligned_query.count('-') + aligned_pattern.count('-')
                    
        
                    if mismatches <= max_mismatches:
                        matches.append({
                            'start': start_position,
                            'end': end_position,
                            'primer': search_primer,
                            'score': alignment.score,
                            'mismatches': mismatches,
                            'gaps': gaps,
                            'sequence': read_name,
                            'is_reverse': search_primer == rev_comp_primer
                        })
                
        return sorted(matches, key=lambda x: (x['start'], -x['score'], x['mismatches']))
    
 
    def find_best_primer_pairs(self, sequence: str, read_name: str) -> List[dict]:

        """
        Find all valid primer pairs in the sequence.
        
        Args:
            sequence: Input DNA sequence
            
        Returns:
            List of dictionaries containing paired sequence information
        """
        #log_message(f"Processing Read={read_name}")
        
        forward_matches = self.smith_waterman_search4(sequence, read_name, self.forward_primer)
        reverse_matches = self.smith_waterman_search4(sequence, read_name, self.reverse_primer)
        
        paired_sequences = []
        used_positions = set()  # Track used positions to avoid overlaps
        
        # Consider all matches as potential primers in either direction
        all_matches = [(pos, 'forward') for pos in forward_matches] + [(pos, 'reverse') for pos in reverse_matches]
        all_matches.sort(key=lambda x: (x[0]['start'], -x[0]['score'], x[0]['mismatches'])) # Sort by position and score and mismatches
  
        
        filtered_matches = []  # List to hold the non-overlapping, filtered matches
        
        for match_data, match_type in all_matches:
            # Extract match start and end positions
            start, end = match_data['start'], match_data['end']
            
            # If the filtered list is empty or there's no overlap, just add the match
            if not filtered_matches:
                filtered_matches.append((match_data, match_type))
            else:
                # Check overlap with the last added match in filtered_matches
                last_match_data, last_match_type = filtered_matches[-1]
                last_start, last_end = last_match_data['start'], last_match_data['end']
                
                # If there is an overlap (last match ends after the current match starts)
                if last_end >= start:
                    # Compare the current match with the last one based on score and mismatches
                    if (match_data['score'], -match_data['mismatches']) > (last_match_data['score'], -last_match_data['mismatches']):
                        # If current match is better, replace the last match
                        filtered_matches[-1] = (match_data, match_type)
                else:
                    # If no overlap, add the current match
                    filtered_matches.append((match_data, match_type))
                    
        filtered_matches.sort(key=lambda x: x[0]['start'])  # Sort by start position
        all_matches = filtered_matches
        log_message(f"Read={read_name} - found {len(all_matches)} primer hits", logging.DEBUG)     
                
        for i, (match1, type1) in enumerate(all_matches[:-1]):  # Look at all pairs of adjacent matches
            
            pos1_start = match1['start']
            pos1_end = match1['end']
            
            # Skip if this position has been used
            if pos1_start in used_positions or pos1_end in used_positions:
                continue
                
            # Look at all subsequent matches
            match2, type2 = all_matches[i+1]
            pos2_start = match2['start']
            pos2_end = match2['end']
            
            # Skip if this position has been used
            if pos2_start in used_positions or pos2_end in used_positions:
                continue
            
            # Set minimum and maximum segment size to consider valid
            distance = pos2_end - pos1_start
            passes_length_quality = 50 <= distance <= 50000
            
             # Check if match1 and match2 form a valid primer pair
            valid_pair = False
            reversed = False
            if type1 == 'forward' and type2 == 'reverse' and match1['is_reverse'] == False and passes_length_quality:
                valid_pair = True
                trimmed_seq = sequence[pos1_start:pos2_end]
            elif type1 == 'reverse' and type2 == 'forward' and match1['is_reverse'] == True and passes_length_quality:
                valid_pair = True
                reversed = True
                match1, match2 = match2, match1
                pos1_start, pos2_start = pos2_start, pos1_start
                pos1_end, pos2_end = pos2_end, pos1_end
                trimmed_seq = str(Seq(sequence[pos2_start:pos1_end]).reverse_complement())
            
            if valid_pair:
                # Valid pair, add both positions to used_positions
                used_positions.update(set(range(pos1_start, pos1_end + 1)) | set(range(pos2_start, pos2_end + 1)))
                # Log the match and trimming information
                log_message(f"Read={read_name} - VALID Primer Pair: ({pos1_start}, {pos1_end}) - ({pos2_start}, {pos2_end}) Length={distance}bp, Reversed: {reversed}")
                
                segment = {
                    'trimmed_seq': trimmed_seq,
                    'forward': (pos1_start, pos1_end),
                    'reverse': (pos2_start, pos2_end),
                    'valid': True,
                    'reversed':reversed,
                    'total_mismatches': match1['mismatches'] + match2['mismatches'],
                    'forward_seq': match1['sequence'],
                    'reverse_seq': match2['sequence']
                }
                 
                paired_sequences.append(segment)
                # Move to the next match for the outer loop to check (i.e., match3)
                continue
            else:
                # Invalid pair, only mark match1 as used
                used_positions.update(set(range(pos1_start, pos1_end + 1)))
                
                # Grab invalid segment for low qual output
                trimmed_seq = sequence[pos1_start:pos2_end]
                log_message(f"Read={read_name} - INVALID Primer Pair: ({pos1_start}, {pos1_end}) - ({pos2_start}, {pos2_end}) Length={distance}bp", logging.DEBUG)
    
                invalid_segment = {
                    'trimmed_seq': trimmed_seq,
                    'forward': (pos1_start, pos1_end),
                    'reverse': (pos2_start, pos2_end),
                    'valid': False,
                    'reversed':reversed,
                    'total_mismatches': match1['mismatches'] + match2['mismatches'],
                    'forward_seq': match1['sequence'],
                    'reverse_seq': match2['sequence']
                }
                paired_sequences.append(invalid_segment)  # Add to the same list with the 'valid' key indicating invalid pair
                
              # Proceed to the next match (i.e., match2 will now be considered with match3)
  
        # Handle the last match if it wasn't paired (no match after it)
        # Change this to if len(all_matches) > 0: to pull any fragments from unpaired primers into lowqual bin (see process records)
        if len(all_matches) == 1:
            last_match, last_type = all_matches[-1]
            last_pos_start = last_match['start']
            last_pos_end = last_match['end']
    
            # If the last match has not been used yet, treat it as an invalid segment
            if last_pos_start not in used_positions and last_pos_end not in used_positions:
                trimmed_seq = sequence[last_pos_start:len(sequence)-1]  # Grab from last match to the end of the sequence
                log_message(f"Read={read_name} - UNPAIRED Primer End Segment: ({last_pos_start}, {last_pos_end}) - Unpaired Primer, end of read, Length={len(trimmed_seq)}bp", logging.DEBUG)
    
                last_segment = {
                    'trimmed_seq': trimmed_seq,
                    'forward': (last_pos_start, last_pos_end),
                    #'reverse': (last_pos_start, last_pos_end),
                    'reverse': (len(sequence)-1, len(sequence)-1),
                    'valid': False,
                    'reversed': False,
                    'total_mismatches': last_match['mismatches'],
                    'forward_seq': last_match['sequence'],
                    'reverse_seq': "NA"
                }
    
                paired_sequences.append(last_segment)

        # If no valid pairs were found, print that the read is binned
        # if not paired_sequences:
        #     log_message(f"Read={read_name} - Binned (No valid primer pairs)", logging.WARNING)
        
        return paired_sequences
    
    def process_record(self, records: List[SeqRecord]) -> List[Tuple[SeqRecord, str]]:
        results_records = []
    
        for record in records:
            seq = str(record.seq)
            read_name = record.id  # Get read name
            log_message(f"Processing read {read_name}", logging.DEBUG)
            matches = self.find_best_primer_pairs(seq, read_name)  # Pass read name
            
            
            #deal with single primer hits - comment this section out if single primer hit segments should be put in lowqual instead
            if len(matches) == 1 and matches[0]['reverse_seq'] == "NA":
               #log_message(f"Read={read_name} - Check - {matches[0]['forward_seq']} Length {len(matches[0]['trimmed_seq'])} Forward {matches[0]['forward']} Reverse {matches[0]['reverse']} Reversed {matches[0]['reversed']}", logging.DEBUG)
               log_message(f"Read={read_name} - Binned (Only one valid primer found). Length {len(matches[0]['trimmed_seq'])} Primer position: {matches[0]['forward']}", logging.DEBUG)
               results_records.append((record, 'binned_1'))
               continue  # Skip to next record
              
            #deal with 0 primer hits
            if not matches:
                # If no matches are found, classify as binned, but always return the record.
                results_records.append((record, 'binned_0'))
                log_message(f"Read={read_name} - Binned (No valid primers found)", logging.DEBUG)
                continue  # Skip to next record
            
            log_message(f"Sort sorting {len(matches)} segments by valid primers and polyA filters.")
            
            for i, match in enumerate(matches, 1):
                
                correct_length = len(match['trimmed_seq']) >= 50 and len(match['trimmed_seq']) <= 50000
                umi_seq = True
                polyT = regex.findall("(TTTTTTTT){e<=0}", str(match['trimmed_seq'])[:150])
                polyT = False
                polyA = regex.findall("(AAAAAAAA){e<=0}", str(match['trimmed_seq'])[-200:])
                valid = match['valid']
                classification = "full_length" if umi_seq and polyA and not polyT and valid and correct_length else "lowqual"
    
                try:
                    ## modify the quality score extraction to allow reversing for flipped reads
                    if match['reversed']:  # If the primer pair is reversed
                        trimmed_qual = record.letter_annotations["phred_quality"][match['reverse'][0]:match['forward'][1]][::-1]
                    else:
                        trimmed_qual = record.letter_annotations["phred_quality"][match['forward'][0]:match['reverse'][1]]
    
                    new_id = f"{record.id}_segment_{i}"
                    new_record = SeqRecord(
                        Seq(match['trimmed_seq']),
                        id=new_id,
                        description=f"{record.description} length={len(match['trimmed_seq'])}",
                        letter_annotations={"phred_quality": trimmed_qual}
                    )
    
                    results_records.append((new_record, classification))
                    log_message(f"Segment={record.id}_segment_{i} - New segment extracted, {classification}", logging.DEBUG)
    
                except Exception:
                    log_message(f"Segment={record.id}_segment_{i} - Skipped segment, could not find quality scores. Slice from {match['forward'][0]} to {match['reverse'][1]}, Length length={len(match['trimmed_seq'])}bp", logging.WARNING)
                    continue
    
        # if not results_records:
        #     print("NO PRIMERS")
        #     # Ensure the 'binned' read is always included if no matches
        #     results_records.extend([(record, 'binned') for record in records])
    
        return results_records

    def worker(self, args: str):
        records_chunk = args
        result = self.process_record(records_chunk)
        
        
        # Count 'binned' directly, and calculate 'classified' by subtracting from total length
        binned_count = sum(1 for _, classification in result if classification == 'binned')
        classified_count = len(result) - binned_count  # Total minus binned

    
        log_message(f"Chunk processed: {len(records_chunk)} reads -> {classified_count} processed segments and {binned_count} binned reads")
    
        return result
        

    def parallel_split_reads(self, input_file: str, processed_output: str, lowqual_output: str, bin_output: str, stats_output: str, num_workers: int, chunk_size: int, verbose: bool):
        """
        Split FastQ reads based on primer pairs and indexes.
        
        Args:
            input_file: Input FastQ file path
            processed_output: Output path for processed reads with UMIs
            lowqual_output: Output path for low quality reads without UMIs and without polyA
            bin_output: Output path for binned reads
            stats_output: Output path for statistics
        """
        stats = {
            'total_sequences': 0,
            'total_segments': 0,
            'processed': 0,
            'lowqual': 0,
            'binned': 0
        }
        stats = {
            'total_reads': 0,
            'processed_reads': 0,
            'total_segments': 0,
            'full_length_segments': 0,
            'lowqual_segments': 0,
            'binned_reads_0': 0,
            'binned_reads_1': 0
        }
        
        processed_records = []
        lowqual_records = []
        binned_records = []
        
        # Open input records                
        if input_file.endswith('.gz'):
            handle = gzip.open(input_file, 'rt')
        else:
            handle = open(input_file, 'r')
        
        with handle:
            records = list(SeqIO.parse(handle, 'fastq'))
            total_records = len(records)
            stats['total_reads'] = total_records

            # Split records into chunks
            #chunk_size = len(records) // num_workers
            chunks = [records[i:i + chunk_size] for i in range(0, len(records), chunk_size)]
            num_workers = min(num_workers, len(chunks))
            
            # Debug print: Total reads and number of chunks
            log_message(f"Processing {total_records} reads in {len(chunks)} chunks.")

            # Prepare arguments for worker processes
            args = [(chunk) for chunk in chunks]

            # Initialize the progress bar
            with Pool(num_workers) as pool:
                with tqdm(total=total_records, desc="Processing Reads", unit=" reads") as pbar:
                    results = []
                    
                    
                    for records_chunk in chunks:  # Iterate over the chunks of records
                        # Pass the chunk to the worker and get the result
                        result = self.worker(records_chunk)
                        results.extend(result)  # Flatten the results
                        
                        # Update the progress bar based on the number of reads in the chunk
                        pbar.update(len(records_chunk))  # Update progress by the size of the chunk

                    # for result in pool.imap(self.worker, args):
                    #     # Flatten the results from the worker
                    #     results.extend(result)
                    #     pbar.update(len(result))  # Update progress bar after processing each chunk

                    all_results = results
                    #stats['total_segments'] = len(all_results)
                    
                    # Debug print: After all chunks are processed
                    log_message(f"Final sorting {stats['total_segments']} segments by valid primers and polyA filters.")
      
                     # Classify and categorize the results
                    for new_record, classification in all_results:
                        if classification == 'full_length':
                            processed_records.append(new_record)
                            stats['full_length_segments'] += 1
                        elif classification == 'lowqual':
                            lowqual_records.append(new_record)
                            stats['lowqual_segments'] += 1
                        elif classification == 'binned_1':
                            binned_records.append(new_record)
                            stats['binned_reads_1'] += 1
                        else:
                            binned_records.append(new_record)
                            stats['binned_reads_0'] += 1
                    # Debug print: After all chunks are processed
                    log_message(f"Sorting segments complete. Writing outputs.")
                    #print(f"\nAll {len(chunks)} chunks processed. Writing output files.")

        stats['processed_reads'] = stats['total_reads'] - (stats['binned_reads_0'] + stats['binned_reads_1']) 
        stats['total_segments'] = stats['full_length_segments'] +  stats['lowqual_segments']
        
        # Write output files
        with gzip.open(processed_output, 'wt') as handle:
            SeqIO.write(processed_records, handle, 'fastq')
            
        with gzip.open(lowqual_output, 'wt') as handle:
            SeqIO.write(lowqual_records, handle, 'fastq')
            
        with gzip.open(bin_output, 'wt') as handle:
            SeqIO.write(binned_records, handle, 'fastq')

        # Write statistics
        with open(stats_output, 'w', newline='') as handle:
            writer = csv.writer(handle)
            writer.writerow(["Metric", "Value"])
            for metric, value in stats.items():
                writer.writerow([metric, value])

        


def setup_logging(results_dir, input_file, verbose):
    """Set up logging based on input file name and verbosity."""
    input_file = os.path.basename(input_file).replace('.fastq.gz', '')
    log_filename = f"output.{input_file}.log"
    log_filename = os.path.join(results_dir, log_filename)
    logging.basicConfig(
        filename=log_filename, 
        level=logging.DEBUG if verbose else logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )

def log_message(message, level=logging.INFO):
    """Log a message if verbose mode is enabled."""
    logging.log(level, message)



        
def main():
    
    parser = argparse.ArgumentParser(description="Split FASTQ files based on index and primer sequences.")
    parser.add_argument("-i", "--input_file", help="Input FASTQ file (gzipped)")
    parser.add_argument("-res", "--results-dir", default=None,
                       help="Directory to store output files (default is basename of input file)")

    parser.add_argument("--processed-output", required=True,
                       help="Output file for processed reads")
    parser.add_argument("--lowqual-output", required=True,
                       help="Output file for low quality reads")
    parser.add_argument("--bin-output", required=True,
                       help="Output file for binned reads")
    parser.add_argument("--stats-output", required=True,
                       help="Output statistics file")
    parser.add_argument("-fp", "--forward-primer",
                       default="AAGCAGTGGTATCAACGCAGAGTGAAT",
                       help="Forward primer sequence")
    parser.add_argument("-rp", "--reverse-primer",
                       default="GTACTCTGCGTTGATACCACTGCTT",
                       help="Reverse primer sequence")
    parser.add_argument("-sw", "--smith-waterman", type=int, choices=[1, 2, 3, 4], default=4,
                        help="Select Smith-Waterman search version (1, 2, 3, or 4). Default is 4.")
    parser.add_argument("-e","--error", type=float, default=0.3, 
                       help="Number of errors allowed, default is 0.3")
    parser.add_argument("--chunk-size", type=int, default=1000,  
                       help="Number of reads per chunk, defualt is 1000")
    parser.add_argument("--num_workers", type=int, default=4,
                       help="Number of parallel workers (CPUs), default is 4")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable detailed logging")

    args = parser.parse_args()
    
    # Set the results directory
    if args.results_dir is None:
        # Default to basename of input file without extension
        args.results_dir = os.path.basename(args.input_file).replace('.fastq.gz', '')
    
    # Make the results directory if it doesn't exist
    if not os.path.exists(args.results_dir):
        os.makedirs(args.results_dir, exist_ok=True)
    
    args.processed_output = os.path.join(args.results_dir, args.processed_output)
    args.lowqual_output = os.path.join(args.results_dir, args.lowqual_output)
    args.bin_output = os.path.join(args.results_dir, args.bin_output)
    args.stats_output = os.path.join(args.results_dir, args.stats_output)
    
    # # Start timer and setup logging using our helper function
    setup_logging(args.results_dir, args.input_file, args.verbose)
    log_message("Script started.")
  
    start_time = time.time()  
    log_message(f"Starting read processing with {args.num_workers} threads and chunk size {args.chunk_size}")
    log_message(f"Scanning for FP {args.forward_primer} and RP {args.reverse_primer} in {args.input_file} with error tolerance {args.error}.")
    splitter = FastqSplitter(args.forward_primer, args.reverse_primer, args.error, args.smith_waterman)
    splitter.parallel_split_reads(
        args.input_file,
        args.processed_output,
        args.lowqual_output,
        args.bin_output,
        args.stats_output,
        args.num_workers,
        args.chunk_size,
        args.verbose
    )
    
    # Calculate and print the total time taken
    end_time = time.time()
    total_time = end_time - start_time

    log_message(f"Outputs written to {args.results_dir}/")
    log_message("Script completed.")
    log_message(f"Total execution time: {total_time:.2f} seconds.")


if __name__ == "__main__":
    main()
