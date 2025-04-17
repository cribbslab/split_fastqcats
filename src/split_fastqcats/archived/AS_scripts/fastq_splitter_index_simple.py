#!/usr/bin/env python3

import sys
import argparse
import gzip
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import csv
from typing import Dict, List, Tuple
import time
from multiprocessing import Pool
from tqdm import tqdm

class FastqSplitter:
    def __init__(self, forward_primer: str, index_dict: Dict[str, str], errors: int):
        
        print("Finding best index matches and splitting segments... \n")
        
        self.forward_primer = forward_primer
        self.index_dict = index_dict
        
        self.max_mismatches = errors
        self.patterns = [self.index_dict[key] + self.forward_primer for key in self.index_dict]
        

    def simple_match_search(self, sequence: str, record_id: str) -> List[dict]:
        
        best_matches = {}

        for pattern in self.patterns:
            pattern_length = len(pattern)
            for i in range(len(sequence) - pattern_length + 1):
                subsequence = sequence[i:i + pattern_length]
                mismatches = sum(1 for a, b in zip(subsequence, pattern) if a != b)

                if mismatches <= self.max_mismatches:
                    if i not in best_matches or mismatches < best_matches[i]['mismatches']:
                        best_matches[i] = {
                            'start': i + (len(subsequence) - pattern_length),
                            'end': i + pattern_length,
                            'mismatches': mismatches,
                            'sequence': subsequence,
                            'pattern': pattern,
                            'record_id': record_id
                        }

        sorted_matches = sorted(best_matches.values(), key=lambda x: x['start'])

        if sorted_matches:
            print(f"Read Name: {record_id}")
            for match in sorted_matches:
                print(f"  Pattern: {match['pattern']}, Start Position: {match['start']}, Mismatches: {match['mismatches']}")

        return sorted_matches
    

    def process_record(self, record: SeqRecord) -> List[Tuple[SeqRecord, str]]:
        seq = str(record.seq)
        matches = self.simple_match_search(seq, record.id)
        
    
        if not matches:
            return [(record, 'binned')]  # No matches found, classify as "binned"
    
        processed_records = []
    
        for i, current_match in enumerate(matches):
            start_position = current_match['start']
            end_position = current_match['end']
            next_start_position = matches[i + 1]['start'] if i + 1 < len(matches) else len(seq)
    
            subsequence = seq[start_position:next_start_position]
    
            # Make sure the subsequence length is within the valid range
            if len(subsequence) < 50 or len(subsequence) > 5000:
                
                continue
    
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
            processed_records.append((new_record, index_label))
    
        return processed_records if processed_records else [(record, 'binned')]  # Classify as "binned" if no processed segments


    def worker(self, records: List[SeqRecord]) -> List[Tuple[SeqRecord, str]]:
        """
        Worker function to process a chunk of FASTQ records.
        """
        results = []
        for record in records:
            processed_segments = self.process_record(record)
            results.extend(processed_segments)
        return results
    
    def parallel_split_reads(self, input_file: str, processed_output: str, bin_output: str, stats_output: str, num_workers: int, chunk_size: int, errors: int):
        stats = {
            'total_sequences': 0,
            'total_segments': 0,
            'processed': 0,
            'binned': 0
        }
    
        processed_records = {str(i): [] for i in range(1, len(self.index_dict))}  # Dictionary to store 12 files, one for each index
        binned_records = []
        index_counts = {sequence: 0 for sequence in self.index_dict.values()}  # Use the actual index sequence as keys
    
    
        start_time = time.time()  # Start timer
        
        with gzip.open(input_file, 'rt') as handle:
            records = list(SeqIO.parse(handle, 'fastq'))
            total_records = len(records)
            stats['total_sequences'] = total_records
    
            # Split records into chunks
            #chunk_size = len(records) // num_workers
            chunks = [records[i:i + chunk_size] for i in range(0, len(records), chunk_size)]
    
            # Prepare arguments for worker processes
            #args = [(chunk) for chunk in chunks]
    
            # Initialize the progress bar
            with Pool(num_workers) as pool:
                with tqdm(total=total_records, desc="Processing Reads", unit=" reads") as pbar:
                    results = []
    
                    for result in pool.imap(self.worker, chunks):
                        # Flatten the results from the worker
                        results.extend(result)
                        pbar.update(len(result))  # Update progress bar after processing each chunk
    
                    all_results = results
                    stats['total_segments'] = len(all_results)
                    print("\n Sorting segments to output files... \n")
    
                    # Classify and categorize the results
                    for new_record, classification in all_results:
                        index_label =  classification
                        if classification in self.index_dict.values():  
                            print(f"Adding to processed: {new_record.id}, Index label: {index_label}, Length {len(new_record)}bp")  # Debugging print
                           
                            # Initialize the list for this index if it doesn't exist
                            if index_label not in processed_records:
                                processed_records[index_label] = []
                            
                            # Add the record to the appropriate index file
                            processed_records[index_label].append(new_record)
                            
                            index_counts[index_label] += 1  # Increment the count for this index
                            stats['processed'] += 1
                        else:
                            binned_records.append(new_record)
                            stats['binned'] += 1



        # Write output files for each index
        for index, records in processed_records.items():
            if records:  # If there are records for this index, write to its file
                output_file = f"{processed_output}_index_{index}.fastq.gz"
                with gzip.open(output_file, 'wt') as handle:
                    SeqIO.write(records, handle, 'fastq')

        # Write binned records
        with gzip.open(bin_output, 'wt') as handle:
            SeqIO.write(binned_records, handle, 'fastq')
    
        # Write statistics
        with open(stats_output, 'w', newline='') as handle:
            writer = csv.writer(handle)
            writer.writerow(["Metric", "Value"])
            for metric, value in stats.items():
                writer.writerow([metric, value])
            # Write the count of segments for each index
            writer.writerow(["Index", "Count"])
            for index, count in index_counts.items():
                writer.writerow([index, count])
    
        # Calculate and print the total time taken
        end_time = time.time()
        total_time = end_time - start_time
        print(f"\nTotal execution time: {total_time:.2f} seconds")
      
      

def parse_indexes(index_args: List[str]) -> Dict[str, str]:
    """Parse index arguments into a dictionary with automatic numbering."""
    index_dict = {}
    for i, sequence in enumerate(index_args, 1):  # Automatically generate label starting from 1
        index_dict[str(i)] = sequence  # Use the sequence as the key and auto-generated label as value
    return index_dict

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
    parser.add_argument("input_file", help="Input FASTQ file (gzipped)")
    parser.add_argument("--processed-output", required=True, help="Output file for processed reads")
    parser.add_argument("--bin-output", required=True, help="Output file for binned reads")
    parser.add_argument("--stats-output", required=True, help="Output statistics file")
    parser.add_argument("-fp", "--forward-primer", default="AAGCAGTGGTATC", help="Forward primer sequence")
    
    # Accept a list of index sequences
    parser.add_argument("--indexes", nargs='+', help="List of index sequences (e.g. 'AAATTTGGGCCC' 'TTTCCCAAAGGG')")
    parser.add_argument("--chunk_size", type=int, default=1000, help="Number of reads per chunk")
    parser.add_argument("--num_workers", type=int, default=4, help="Number of parallel workers (CPUs), default is 4")
    parser.add_argument("--errors", type=int, default=4, help="Number of errors allowed, default is 4")

    args = parser.parse_args()

    # Use the new parse_indexes function
    index_dict = parse_indexes(args.indexes) if args.indexes else default_indexes
    
    splitter = FastqSplitter(args.forward_primer, index_dict, args.errors)
    splitter.parallel_split_reads(
        args.input_file,
        args.processed_output,
        args.bin_output,
        args.stats_output,
        args.num_workers,
        args.chunk_size,
        args.errors
    )
    
if __name__ == "__main__":
    main()

