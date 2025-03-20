import sys
import argparse
import gzip
import regex
import math
import itertools
import parasail  # Fast SIMD-based sequence alignment
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from concurrent.futures import ProcessPoolExecutor
from tqdm import tqdm
from typing import Dict, List, Tuple, Union

class FastqSplitter:
    def __init__(self, forward_primer: str, index_dict: Dict[str, str], errors: int):
        print("Finding best index matches and splitting segments...\n")
        self.forward_primer = forward_primer
        self.index_dict = index_dict
        self.max_mismatches = errors
        self.patterns = [self.index_dict[key] + self.forward_primer for key in self.index_dict]
        self.match_score = 2
        self.mismatch_score = -1
        self.open_gap_score = 3
        self.extend_gap_score = 1

    # def smith_waterman_search(self, sequence: str, record_id: str, errors: Union[int, float] = None) -> List[dict]:
    #     """
    #     Perform Smith-Waterman local alignment to find the best pattern (index + primer) match in the sequence.
    # 
    #     Args:
    #         sequence: The input DNA sequence.
    #         record_id: The identifier of the read.
    #         errors: If float (0 < errors < 1), it's used as a similarity threshold factor.
    #                 If integer, it's used as the max mismatches allowed.
    # 
    #     Returns:
    #         A sorted list of match dictionaries with alignment info.
    #     """
    #     best_matches = []
    #     if errors is None:
    #         errors = self.max_mismatches  # Use the user-supplied value from __init__
    # 
    # 
    #     for pattern in self.patterns:
    # 
    #         pattern_length = len(pattern)
    # 
    #         # If `errors` is a float (0 < errors < 1), it's treated as a similarity threshold factor
    #         if isinstance(errors, float) and errors < 1:
    #             min_score_threshold = (1 - errors) * pattern_length * self.match_score
    #             max_mismatches = math.floor(errors * pattern_length)  # No strict mismatch limit
    #         else:  # Treat `errors` as an integer for the maximum mismatches
    #             errors = round(errors)  # Round the float to an integer
    #             min_score_threshold = (pattern_length - errors) * self.match_score
    #             max_mismatches = errors
    # 
    # 
    #         #alignment = parasail.sw_trace(sequence, pattern, 2, 1, parasail.blosum62)
    #         alignment = parasail.sw_trace_striped_16(
    #             sequence, pattern,
    #             #self.match_score,   # Match score
    #             #self.mismatch_score, # Mismatch penalty
    #             self.open_gap_score, # Gap opening penalty
    #             self.extend_gap_score, # Gap extension penalty
    #             parasail.matrix_create("ACGT", self.match_score, self.mismatch_score)
    #             #parasail.dnafull  # DNA-specific scoring matrix
    #         )
    #         score = alignment.score
    #         mismatches = len(pattern) - (score // 2)  # Approximate mismatches
    #         #if score >= (len(pattern) - self.max_mismatches) * 2:
    #         if score >= min_score_threshold:
    #             start_position = alignment.end_query - len(pattern) + 1
    #             # Log the start position, score, and mismatches
    #             #print(f"DEBUG: Record={record_id}, Pattern={pattern}, Start={start_position}, Mismatches={mismatches}, Score={score}")
    #             best_matches.append({
    #                 'start': start_position,
    #                 'score': score,
    #                 'pattern': pattern,
    #                 'mismatches': mismatches,
    #                 'record_id': record_id,
    #                 'min_score_threshold': min_score_threshold
    #             })
    #         
    #     return sorted(best_matches, key=lambda x: (x['start'], -x['score']))

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
        best_matches = []
        if errors is None:
            errors = self.max_mismatches  # Default to class-level setting

        for pattern in self.patterns:
            pattern_length = len(pattern)

            # If `errors` is a float (0 < errors < 1), it's treated as a similarity threshold factor
            if isinstance(errors, float) and errors < 1:
                min_score_threshold = (1 - errors) * pattern_length * self.match_score
                max_mismatches = math.floor(errors * pattern_length)  # No strict mismatch limit
            else:  # Treat `errors` as an integer for the maximum mismatches
                errors = round(errors)  # Round the float to an integer
                min_score_threshold = (pattern_length - errors) * self.match_score
                max_mismatches = errors

            # Perform Parasail alignment
            alignment = parasail.sw_trace_striped_16(
                sequence, pattern,
                self.open_gap_score, # Gap opening penalty
                self.extend_gap_score, # Gap extension penalty
                parasail.matrix_create("ACGT", self.match_score, self.mismatch_score)
                #parasail.dnafull  # DNA-specific scoring matrix
            )

            # **Check if alignment score is valid**
            if alignment is None or alignment.score is None:
                continue  # Skip this pattern if alignment failed

            # Extract actual alignment details
            aligned_query = alignment.traceback.query
            aligned_pattern = alignment.traceback.ref

            # Count mismatches and gaps
            mismatches = sum(1 for q, p in zip(aligned_query, aligned_pattern) if q != p and q != '-' and p != '-')
            gaps = sum(1 for q, p in zip(aligned_query, aligned_pattern) if q == '-' or p == '-')

            # Ensure start position is valid
            start_position = max(0, alignment.end_query - pattern_length + 1)

            # Allow a maximum of `max_mismatches` mismatches (including gaps)
            if mismatches + gaps <= max_mismatches and alignment.score >= min_score_threshold:
                #print(f"DEBUG: Record={record_id}, Pattern={pattern}, Start={start_position}, Mismatches={mismatches}, Score={alignment.score}")

                best_matches.append({
                    'start': start_position,
                    'score': alignment.score,
                    'pattern': pattern,
                    'mismatches': mismatches,
                    'gaps': gaps,
                    'record_id': record_id,
                    'min_score_threshold': min_score_threshold
                })

        # Sort by best match (highest score, lowest mismatches, then lowest start position)
        return sorted(best_matches, key=lambda x: (x['start'], x['mismatches'], x['gaps'], -x['score']))

    def process_record(self, record: SeqRecord) -> List[Tuple[SeqRecord, str]]:
        seq = str(record.seq)
        matches = self.smith_waterman_search(seq, record.id)
        if not matches:
            print(f"DEBUG: Read={record.id} -> Binned")
            return [(record, 'binned')]
        
        segments = []
        for i, match in enumerate(matches):
        #for match in matches:
            start_position = match['start']
            next_start_position = matches[i + 1]['start'] if i + 1 < len(matches) else len(seq)
            

              
            subsequence = seq[start_position:next_start_position]
            
            # Make sure the subsequence length is within the valid range
            if len(subsequence) < 50 or len(subsequence) > 10000:
                print(f"DEBUG: Skipping segment {record.id}_segment_{i+1}, Length: {len(subsequence)}")
                continue
            
            index_label = match['pattern'].split(self.forward_primer)[0]
            new_id = f"{record.id}_segment_{len(segments) + 1}"
            trimmed_qual = record.letter_annotations.get("phred_quality", [])
            trimmed_qual = trimmed_qual[start_position:next_start_position] if trimmed_qual else []
    
            #trimmed_qual = record.letter_annotations.get("phred_quality", [])[start_position:next_start_position] if record.letter_annotations else []
            new_record = SeqRecord(
                Seq(subsequence), id=new_id, description=f"{record.description} length={len(subsequence)}",
                letter_annotations={"phred_quality": trimmed_qual} if trimmed_qual else {}
            )
            
            print(f"DEBUG: Read={record.id}, Segment={new_id} -> Assigned to {index_label}, Start={start_position}, Gaps={match['gaps']}, Mismatches={match['mismatches']}, Score={match['score']}, Score threshold: {match['min_score_threshold']}, Length={len(subsequence)}bp")
            #print(f"DEBUG: Segments returned for {record.id}: {[(seg[0].id, seg[1]) for seg in segments]}")

            segments.append((new_record, index_label))

        return segments

    def worker(self, records: List[SeqRecord]) -> List[Tuple[SeqRecord, str]]:
        return [segment for record in records for segment in self.process_record(record)]

    def process_file_in_chunks(self, input_file: str, chunk_size: int):
        with gzip.open(input_file, 'rt') as handle:
            records = list(SeqIO.parse(handle, 'fastq'))
            print(f"Total Reads: {len(records)}")  # Print total reads for visibility
            print(f"Processing in: {chunk_size} read chunks\n")  # Print total reads for visibility
            for i in range(0, len(records), chunk_size):
                yield records[i:i + chunk_size]

    def parallel_split_reads(self, input_file: str, processed_output: str, bin_output: str, stats_output: str, num_workers: int, chunk_size: int):
        stats = {'total_reads': 0, 'processed_reads': 0, 'total_segments': 0, 'processed_segments': 0, 'binned_reads': 0}
        processed_records = {key: [] for key in self.index_dict.values()}
        binned_records = []
        index_counts = {sequence: 0 for sequence in self.index_dict.values()}  # Use the actual index sequence as keys
        
        

        with ProcessPoolExecutor(max_workers=num_workers) as executor:
            with tqdm(total=stats['total_reads'], desc="Processing Reads", unit=" reads",  dynamic_ncols=True) as pbar:
                for chunk in self.process_file_in_chunks(input_file, chunk_size):
                    stats['total_reads'] += len(chunk)
                    results = list(executor.map(self.worker, [chunk]))
        
                    for res in results:
                        for new_record, classification in res:
                            stats['total_segments'] += 1
                            if classification in processed_records:
                                processed_records[classification].append(new_record)
                                index_counts[classification] += 1
                                stats['processed_segments'] += 1
                            else:
                                binned_records.append(new_record)
                                stats['binned_reads'] += 1
                    
                    # Update progress bar by number of reads processed, not chunks
                    pbar.update(len(chunk))
                    # Count processed reads correctly (total - binned)
                    stats['processed_reads'] = stats['total_reads'] - stats['binned_reads']


        # with ProcessPoolExecutor(max_workers=num_workers) as executor:
        #     for chunk in tqdm(self.process_file_in_chunks(input_file, chunk_size), desc="Processing Reads", unit=" reads"):
        #         stats['total_reads'] += len(chunk)
        #         results = list(executor.map(self.worker, [chunk]))
        #         for res in results:
        #             for new_record, classification in res:
        #                 stats['total_segments'] += 1
        #                 if classification in processed_records:
        #                     processed_records[classification].append(new_record)
        #                     index_counts[classification] += 1  # Increment the count for this index
        #                     stats['processed_segments'] += 1
        #                 else:
        #                     binned_records.append(new_record)
        #                     stats['binned_reads'] += 1
        #         # Count processed reads correctly (total - binned)
        #         stats['processed_reads'] = stats['total_reads'] - stats['binned_reads']
                
        for index, records in processed_records.items():
            #print(f"DEBUG: Writing {len(records)} segments to {processed_output}_index_{index}.fastq.gz")
            if records:
              
                with gzip.open(f"{processed_output}_index_{index}.fastq.gz", 'wt') as handle:
                    SeqIO.write(records, handle, 'fastq')

        with gzip.open(bin_output, 'wt') as handle:
            SeqIO.write(binned_records, handle, 'fastq')

        stats['processed_reads'] = stats['total_reads'] - stats['binned_reads']
        
        with open(stats_output, 'w') as handle:
            handle.write("Metric,Value\n")
            for key, value in stats.items():
                handle.write(f"{key},{value}\n")
            # Write the count of segments for each index
            handle.write("\nIndex,SegmentCount\n")
            for index, count in index_counts.items():
                handle.write(f"{index},{count}\n")
        
        print(f"\nCompleted analysing {stats['total_reads']} reads\n")


def parse_indexes(index_args: List[str]) -> Dict[str, str]:
    return {str(i): sequence for i, sequence in enumerate(index_args, 1)}

def main():
    parser = argparse.ArgumentParser(description="Optimized FASTQ Splitter")
    parser.add_argument("input_file", help="Input FASTQ file (gzipped)")
    parser.add_argument("--processed-output", required=True, help="Output file for processed reads")
    parser.add_argument("--bin-output", required=True, help="Output file for binned reads")
    parser.add_argument("--stats-output", required=True, help="Output statistics file")
    parser.add_argument("-fp", "--forward-primer", default="AAGCAGTGGTATC", help="Forward primer sequence")
    parser.add_argument("--indexes", nargs='+', help="List of index sequences")
    parser.add_argument("--num_workers", type=int, default=4, help="Number of parallel workers")
    parser.add_argument("-e","--errors", type=float, default=2, help="Max allowed errors")
    parser.add_argument("--chunk-size", type=int, default=1000, help="Number of reads per chunk")
    args = parser.parse_args()
    
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

    
    index_dict = parse_indexes(args.indexes) if args.indexes else default_indexes
    splitter = FastqSplitter(args.forward_primer, index_dict, args.errors)
    splitter.parallel_split_reads(args.input_file, args.processed_output, args.bin_output, args.stats_output, args.num_workers, args.chunk_size)

if __name__ == "__main__":
    main()
