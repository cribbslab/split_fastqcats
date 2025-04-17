#!/usr/bin/env python3
import os
import logging
import sys
import argparse
import gzip
import math
from Bio import SeqIO, Align
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import parasail  # Fast SIMD-based sequence alignment
import csv
from typing import Dict, List, Tuple, Union
import time
from multiprocessing import Pool
from tqdm import tqdm

class FastqSplitter:
    def __init__(self, forward_primer: str, index_dict: Dict[str, str], mismatches: int):
        
        log_message("Finding best barcode matches and splitting contcatenated reads into segments...")
        
        self.forward_primer = forward_primer[:10]
        self.index_dict = index_dict
        self.max_mismatches = mismatches
        self.patterns = [self.index_dict[key] + self.forward_primer for key in self.index_dict]
        """Configure the Smith-Waterman aligner with scoring parameters."""
        #Parasail arguments
        self.match_score = 2
        self.mismatch_score = -1
        self.open_gap_score = 1
        self.extend_gap_score = 1
    

    #Parasail aligner tool with 100bp window and 50bp overlapping step - very fast -- 280 reads/sec
    def smith_waterman_search(self, sequence: str, record_id: str, errors: Union[int, float] = None) -> List[dict]:
        """
        Perform Smith-Waterman local alignment to find the best pattern (index + primer) match in the sequence.

        Args:
            sequence: The input DNA sequence.
            record_id: The identifier of the read.
            errors: If float (0 < errors < 1), it's used as a similarity threshold factor.
                    If integer, it's used as the max mismatches allowed.

        Returns:
            A sorted list of match dictionaries with alignment info.
        """
        
        #log_message(f"Performing Smith-Waterman search for Read={record_id} with primer={self.forward_primer} and barcodes {list(self.index_dict.values())}", logging.DEBUG)

        best_matches = []
        pattern_length = len(self.patterns[0])
        seen_matches = set()
        
        # Determine min_score_threshold
        errors = self.max_mismatches  # Default to class-level setting
        # If `errors` is a float (0 < errors < 1), it's treated as a similarity threshold factor
        if isinstance(errors, float) and errors < 1:
            min_score_threshold = (1 - errors) * pattern_length * self.match_score
            max_mismatches = math.floor(errors * pattern_length)  # No strict mismatch limit
        else:  # Treat `errors` as an integer for the maximum mismatches
            errors = round(errors)  # Round the float to an integer
            min_score_threshold = (pattern_length - errors) * self.match_score
            max_mismatches = errors
        
        window_size = 100
        step_size = 50
        
        for start in range(0, len(sequence) - window_size + 1, step_size):
            window = sequence[start:min(start + window_size, len(sequence))]  # Extract window
    
    
            for pattern in self.patterns:
                
                # Perform Parasail alignment
                alignment = parasail.sw_trace_striped_16(
                    window, pattern,
                    self.open_gap_score, # Gap opening penalty
                    self.extend_gap_score, # Gap extension penalty
                    parasail.matrix_create("ACGT", self.match_score, self.mismatch_score)
                    #parasail.dnafull  # DNA-specific scoring matrix
                )
                
                # **Check if alignment score is valid**
                if alignment is None or alignment.score is None:
                    print(alignment)
                    continue  # Skip this pattern if alignment failed
    
                
                # Ensure start position is valid
                start_position = max(0, alignment.end_query - pattern_length + 1 + start)
                end_position = min(len(sequence), start_position + pattern_length)
                
                
                # Create a unique identifier for each match
                match_key = (pattern, start_position)
          
                # Avoid adding duplicates
                if match_key in seen_matches:
                    continue  # Skip duplicates
          
                # Otherwise, add to results
                seen_matches.add(match_key)
    
                # Allow a maximum of `max_mismatches` mismatches (including gaps)
                if alignment.score >= min_score_threshold:
                    
                    # Extract actual alignment details
                    aligned_query = alignment.traceback.query
                    aligned_pattern = alignment.traceback.ref
                    
                    # Count mismatches and gaps
                    #counts mismatches and gaps together
                    mismatches = sum(1 for a, b in zip(aligned_query, aligned_pattern) if a != b)
                    
                    #counts mismatch and gaps separately
                    #mismatches = sum(1 for q, p in zip(aligned_query, aligned_pattern) if q != p and q != '-' and p != '-')
                    gaps = aligned_query.count('-') + aligned_pattern.count('-')
                    
                    if mismatches <= max_mismatches:
                        #print(f"DEBUG: Record={record_id}, Pattern={pattern}, Start={start_position}, Mismatches={mismatches}, Score={alignment.score}")
                        log_message(f"Record={record_id}, Pattern={pattern}, Start={start_position+start}, Mismatches={mismatches}, Score={alignment.score}", logging.DEBUG)
    
                        best_matches.append({
                            'start': start_position,
                            'end': end_position,
                            'score': alignment.score,
                            'pattern': pattern,
                            'mismatches': mismatches,
                            #'gaps': gaps, #uncomment if want to count gaps + mismatches separately
                            'record_id': record_id,
                            'min_score_threshold': min_score_threshold
                        })
            
        # Sort by best match criteria (lowest mismatches, highest score, earliest start position)
        sorted_matches = sorted(best_matches, key=lambda x: (x['start'], -x['score'], x['mismatches']))
  
        filtered_matches = []  # List to hold the non-overlapping, filtered matches
        
        for match_data in sorted_matches:
            # Extract match start and end positions
            start, end = match_data['start'], match_data['end']
            
            # Sort by best match (highest score, lowest mismatches, then lowest start position)
            # If the filtered list is empty or there's no overlap, just add the match
            if not filtered_matches:
                filtered_matches.append(match_data)
            else:
                # Check overlap with the last added match in filtered_matches
                last_match_data = filtered_matches[-1]
                last_start, last_end = last_match_data['start'], last_match_data['end']
                
                # If there is an overlap (last match ends after the current match starts)
                if last_end >= start:
                    # Compare the current match with the last one based on score and mismatches
                    if (match_data['score'], -match_data['mismatches']) > (last_match_data['score'], -last_match_data['mismatches']):
                        # If current match is better, replace the last match
                        filtered_matches[-1] = (match_data)
                else:
                    # If no overlap, add the current match
                    filtered_matches.append(match_data)
                    
        filtered_matches.sort(key=lambda x: x['start'])  # Sort by start position
        sorted_matches = filtered_matches
        log_message(f"Read={record_id} - found {len(sorted_matches)} primer hits", logging.DEBUG) 
      
        return sorted_matches



    def process_record(self, records: List[SeqRecord]) -> List[Tuple[SeqRecord, str]]:
        
        results_records = []
        
        for record in records:
            seq = str(record.seq)
            read_name = record.id  # Get read name
            log_message(f"Processing read {read_name}", logging.DEBUG)
            matches = self.smith_waterman_search(seq, read_name, self.max_mismatches)  # Pass read name
    
            if not matches:
                 # If no matches are found, classify as binned, but always return the record.
                log_message(f"Read={record.id} -> {len(matches)} barcode/primer hits found -> Binned", logging.DEBUG)
                results_records.append((record, 'binned'))
                continue  # Skip to next record
        
            processed_records = []
            
            for i, current_match in enumerate(matches):
                start_position = current_match['start']
                end_position = current_match['end']
                next_start_position = matches[i + 1]['start'] if i + 1 < len(matches) else len(seq)
                
                subsequence = seq[start_position:next_start_position]
        
                new_id = f"{record.id}_segment_{i + 1}"
        
                trimmed_qual = record.letter_annotations.get("phred_quality", [])
                trimmed_qual = trimmed_qual[start_position:next_start_position] if trimmed_qual else []
        
                new_record = SeqRecord(
                    Seq(subsequence),
                    id=new_id,
                    description=f"{record.description} length={len(subsequence)}",
                    letter_annotations={"phred_quality": trimmed_qual} if trimmed_qual else {}
                )
        
                # Now, we add the record to the correct index's processed output
                index_label = current_match['pattern'].split(self.forward_primer)[0]  # Get the index label (before the primer)
        
                # Add the record to the processed list with the corresponding index label
                # Make sure the subsequence length is within the valid range
                if len(subsequence) < 50 or len(subsequence) > 50000:
                    log_message(f"Read={record.id} -> {len(matches)} barcode/primer hits found, Skipping segment {new_id}, Length: {len(subsequence)}", logging.DEBUG)
                    results_records.append((new_record, "lowqual"))
                else:
                    results_records.append((new_record, index_label))
                    log_message(f"Read={record.id} -> {len(matches)} barcode/primer hits found, Segment={new_id} -> Assigned to {index_label}, Start={start_position}, Mismatches={current_match['mismatches']}, Score={current_match['score']}, Score threshold: {current_match['min_score_threshold']}, Length={len(subsequence)}bp", logging.DEBUG)
           
    
        return results_records

    def worker(self, args: str):
        """
        Worker function to process a chunk of FASTQ records.
        """
        records_chunk = args
        result = self.process_record(records_chunk)
        # Count 'binned' directly, and calculate 'classified' by subtracting from total length
        binned_count = sum(1 for _, classification in result if classification == 'binned')
        classified_count = len(result) - binned_count  # Total minus binned

        log_message(f"Chunk processed: {len(records_chunk)} reads -> {classified_count} processed segments and {binned_count} binned reads")
    
        return result
    
    def parallel_split_reads(self, input_file: str, processed_output: str, lowqual_output: str, bin_output: str, stats_output: str, num_workers: int, chunk_size: int):
        stats = {
            'total_reads': 0,
            'processed_reads': 0,
            'total_segments': 0,
            'processed_segments': 0,
            'lowqual_segments': 0,
            'binned_reads': 0
        }
    
        processed_records = {str(i): [] for i in range(1, len(self.index_dict))}  # Dictionary to store 12 files, one for each index
        index_counts = {sequence: 0 for sequence in self.index_dict.values()}  # Use the actual index sequence as keys
        lowqual_records = []
        binned_records = []
        # Initialize the list for this index if it doesn't exist
        for index_label in self.index_dict.values():
            if index_label not in processed_records:
                processed_records[index_label] = []
        
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
            chunks = [records[i:i + chunk_size] for i in range(0, len(records), chunk_size)]
            num_workers = min(num_workers, len(chunks))
    
             # Debug print: Total reads and number of chunks
            log_message(f"Processing {total_records} reads in {len(chunks)} chunks.")
            
            # Prepare arguments for worker processes
            args = [(chunk) for chunk in chunks]
    
            # Initialize the progress bar
            with Pool(num_workers) as pool:
                with tqdm(total=total_records, desc="Processing Reads", unit=" reads", dynamic_ncols=True,  unit_scale=True) as pbar:
                    results = []
                    
                    for records_chunk in chunks:  # Iterate over the chunks of records
                        # Pass the chunk to the worker and get the result
                        result = self.worker(records_chunk)
                        results.extend(result)  # Flatten the results
                        
                        # Update the progress bar based on the number of reads in the chunk
                        pbar.update(len(records_chunk))  # Update progress by the size of the chunk

    
                  
                    stats['total_segments'] = len(results)
                    
                     # Debug print: After all chunks are processed
                    log_message(f"Preparing {stats['total_segments']} segments by barcode for writing outputs.")
      
                    # Classify and categorize the results
                    for new_record, classification in results:
                        if classification == 'binned':
                            binned_records.append(new_record)
                            stats['binned_reads'] += 1
                        
                        elif classification in self.index_dict.values():
                           
                            #log_message(f"Adding to processed: {new_record.id}, Index label: {classification}", logging.DEBUG)  # Debugging print
                          
                            # Add the record to the appropriate index file
                            processed_records[classification].append(new_record)
                            index_counts[classification] += 1  # Increment the count for this index
                            stats['processed_segments'] += 1
                        else:
                            lowqual_records.append(new_record)
                            stats['lowqual_segments'] += 1
                            
                    # Debug print: After all chunks are processed
                    log_message(f"Sorting segments complete. Writing outputs.")

        # Write output files for each index
        for index, records in processed_records.items():
            if records:  # If there are records for this index, write to its file
                output_file = f"{processed_output}_index_{index}.fastq.gz"
                with gzip.open(output_file, 'wt') as handle:
                    SeqIO.write(records, handle, 'fastq')

        # Write binned records
        with gzip.open(bin_output, 'wt') as handle:
            SeqIO.write(binned_records, handle, 'fastq')
        
        stats['processed_reads'] = stats['total_reads'] - stats['binned_reads']
        stats['total_segments'] = stats['processed_segments'] +  stats['lowqual_segments']
        # Write lowqual records    
        with gzip.open(lowqual_output, 'wt') as handle:
            SeqIO.write(lowqual_records, handle, 'fastq')
    
        # Write statistics
        with open(stats_output, 'w', newline='') as handle:
            writer = csv.writer(handle)
            writer.writerow(["Metric", "Value"])
            for metric, value in stats.items():
                writer.writerow([metric, value])
            # Write the count of segments for each index
            writer.writerow(["Index", "SegmentCount"])
            for index, count in index_counts.items():
                writer.writerow([index, count])
    
        
      
      

def parse_indexes(index_args: List[str]) -> Dict[str, str]:
    """Parse index arguments into a dictionary with automatic numbering."""
    index_dict = {}
    for i, sequence in enumerate(index_args, 1):  # Automatically generate label starting from 1
        index_dict[str(i)] = sequence  # Use the sequence as the key and auto-generated label as value
    return index_dict


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
    
    default_indexes = {
        '1': 'AAATTTGGGCCC',
        '2': 'TTTCCCAAAGGG',
        '3': 'GGGAAACCCTTT',
        '4': 'CCCGGGTTTAAA',
        '5': 'AAACCCGGGAAA',
        '6': 'TTTGGGAAATTT',
        '7': 'GGGTTTCCCGGG',
        '8': 'CCCAAATTTCCC',
        '9': 'AAAGGGAAAGGG',
        '10': 'TTTAAATTTAAA',
        '11': 'GGGCCCGGGCCC',
        '12': 'CCCTTTCCCTTT'
    }

    
    parser = argparse.ArgumentParser(description="Split FASTQ files based on index and primer sequences.")
    parser.add_argument("-i", "--input_file", help="Input FASTQ file (gzipped)")
    parser.add_argument("-res", "--results-dir", default=None,
                       help="Directory to store output files (default is basename of input file)")

    parser.add_argument("--processed-output", required=True, help="Output file for processed reads")
    parser.add_argument("--lowqual-output", required=True,
                       help="Output file for low quality reads")
    parser.add_argument("--bin-output", required=True, help="Output file for binned reads")
    parser.add_argument("--stats-output", required=True, help="Output statistics file")
    parser.add_argument("-fp", "--forward-primer", default="AAGCAGTGGTATCAACGCAGAGT", help="Forward primer sequence")
    
    # Accept a list of index sequences
    parser.add_argument("--indexes", nargs='+', help="List of index sequences (e.g. 'AAATTTGGGCCC' 'TTTCCCAAAGGG')")
    parser.add_argument("--chunk_size", type=int, default=1000, help="Number of reads per chunk")
    parser.add_argument("--num_workers", type=int, default=4, help="Number of parallel workers (CPUs), default is 4")
    parser.add_argument("-e", "--error", type=float, default=3, help="Number of errors allowed, default is 4")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable detailed logging")
    args = parser.parse_args()

    # Use the new parse_indexes function
    index_dict = parse_indexes(args.indexes) if args.indexes else default_indexes
    
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
    start_time = time.time()  
    log_message("Script started.")
    log_message(f"Starting read processing with {args.num_workers} threads and chunk size {args.chunk_size}")
    log_message(f"Scanning for FP {args.forward_primer} and barcodes {list(index_dict.values())} in {args.input_file} with error tolerance {args.error}.")
    
    splitter = FastqSplitter(args.forward_primer, index_dict, args.error)
    splitter.parallel_split_reads(
        args.input_file,
        args.processed_output,
        args.lowqual_output,
        args.bin_output,
        args.stats_output,
        args.num_workers,
        args.chunk_size
    )
    
    # Calculate and print the total time taken
    end_time = time.time()
    total_time = end_time - start_time

    log_message(f"Outputs written to {args.results_dir}/")
    log_message("Script completed.")
    log_message(f"Total execution time: {total_time:.2f} seconds.")
    
if __name__ == "__main__":
    main()

