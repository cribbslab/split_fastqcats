#!/usr/bin/env python3

"""
FastQ Splitter - A tool for processing and splitting FastQ reads based on primer sequences.

This script processes FastQ files containing concatenated reads and splits them based on
primer sequences, UMIs, and polyA tails. It implements quality control measures and
provides detailed statistics about the processing.
"""

import sys
import argparse
import gzip
from Bio import SeqIO, Align
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import csv
from typing import Dict, List, Tuple, Optional

class FastqSplitter:
    def __init__(self, forward_primer: str, reverse_primer: str, index_dict: Dict[str, List[str]]):
        """
        Initialize the FastQ splitter with primers and index dictionary.
        
        Args:
            forward_primer: Forward primer sequence
            reverse_primer: Reverse primer sequence
            index_dict: Dictionary mapping index numbers to [start_index, end_index] pairs
        """
        self.forward_primer = forward_primer
        self.reverse_primer = reverse_primer
        self.index_dict = index_dict
        self.aligner = Align.PairwiseAligner()
        self._setup_aligner()

    def _setup_aligner(self):
        """Configure the Smith-Waterman aligner with scoring parameters."""
        self.aligner.mode = 'local'
        self.aligner.match_score = 2
        self.aligner.mismatch_score = -1
        self.aligner.open_gap_score = -1
        self.aligner.extend_gap_score = -0.5

    def smith_waterman_search(self, sequence: str, primer: str, max_mismatches: int = 5) -> List[dict]:
        """
        Perform Smith-Waterman local alignment to find primer matches.
        
        Args:
            sequence: Input DNA sequence
            primer: Primer sequence to search for
            max_mismatches: Maximum allowed mismatches
            
        Returns:
            List of dictionaries containing match information
        """
        matches = []
        primer_length = len(primer)
        min_score_threshold = 0.60 * primer_length * self.aligner.match_score  

        # Also search for reverse complement
        rev_comp_primer = str(Seq(primer).reverse_complement())
        
        for search_primer in [primer, rev_comp_primer]:
            for i in range(len(sequence) - primer_length + 1):
                sub_seq = sequence[i:i + primer_length]
                alignments = self.aligner.align(sub_seq, search_primer)
                
                if not alignments:
                    continue
                    
                alignment = alignments[0]
                score = alignment.score
                
                if score >= min_score_threshold:
                    mismatches = sum(1 for a, b in zip(sub_seq, search_primer) if a != b)
                    if mismatches <= max_mismatches:
                        matches.append({
                            'start': i,
                            'end': i + primer_length,
                            'score': score,
                            'mismatches': mismatches,
                            'sequence': sub_seq,
                            'is_reverse': search_primer == rev_comp_primer
                        })
                    
        return sorted(matches, key=lambda x: (x['mismatches'], -x['score'], x['start']))

    def find_best_primer_pairs(self, sequence: str) -> List[dict]:
        """
        Find all valid primer pairs in the sequence.
        
        Args:
            sequence: Input DNA sequence
            
        Returns:
            List of dictionaries containing paired sequence information
        """
        forward_matches = self.smith_waterman_search(sequence, self.forward_primer)
        reverse_matches = self.smith_waterman_search(sequence, self.reverse_primer)
        
        paired_sequences = []
        used_positions = set()  # Track used positions to avoid overlaps
        
        # Consider all matches as potential primers in either direction
        all_matches = [(pos, 'forward') for pos in forward_matches] + [(pos, 'reverse') for pos in reverse_matches]
        all_matches.sort(key=lambda x: x[0]['start'])  # Sort by position
        
        for i, (match1, type1) in enumerate(all_matches[:-1]):  # Look at all pairs of adjacent matches
            pos1_start = match1['start']
            pos1_end = match1['end']
            
            # Skip if this position has been used
            if pos1_start in used_positions or pos1_end in used_positions:
                continue
                
            for match2, type2 in all_matches[i+1:]:  # Look at all subsequent matches
                pos2_start = match2['start']
                pos2_end = match2['end']
                
                # Skip if this position has been used
                if pos2_start in used_positions or pos2_end in used_positions:
                    continue
                
                # Check if we have a valid pair (one forward, one reverse)
                valid_pair = False
                if type1 == 'forward' and type2 == 'reverse':
                    valid_pair = True
                elif type1 == 'reverse' and type2 == 'forward':
                    # Swap them to maintain forward->reverse order
                    match1, match2 = match2, match1
                    pos1_start, pos2_start = pos2_start, pos1_start
                    pos1_end, pos2_end = pos2_end, pos1_end
                    valid_pair = True
                
                if valid_pair:
                    # Check if the distance between primers is reasonable
                    distance = pos2_start - pos1_end
                    if 10 <= distance <= 2000:
                        segment = {
                            'trimmed_seq': sequence[pos1_end:pos2_start],
                            'forward': (pos1_start, pos1_end),
                            'reverse': (pos2_start, pos2_end),
                            'total_mismatches': match1['mismatches'] + match2['mismatches'],
                            'forward_seq': match1['sequence'],
                            'reverse_seq': match2['sequence']
                        }
                        paired_sequences.append(segment)
                        used_positions.update([pos1_start, pos1_end, pos2_start, pos2_end])
                        break  # Move to next first position
        
        return paired_sequences

    def process_record(self, record: SeqRecord) -> List[Tuple[SeqRecord, str]]:
        """
        Process a single FastQ record and split it into multiple records if needed.
        
        Args:
            record: SeqRecord object to process
            
        Returns:
            List of tuples (processed record, classification)
        """
        seq = str(record.seq)
        matches = self.find_best_primer_pairs(seq)
        
        if not matches:
            return [(record, 'binned')]

        processed_records = []
        for i, match in enumerate(matches, 1):
            trimmed_seq = match['trimmed_seq']
            
            # Skip segments that are too short or too long
            if len(trimmed_seq) < 50 or len(trimmed_seq) > 5000:
                continue
                
            # Create new record ID with segment number
            new_id = f"{record.id}_segment_{i}"
            
            try:
                # Get quality scores for the trimmed segment
                trimmed_qual = record.letter_annotations["phred_quality"][
                    match['forward'][1]:match['reverse'][0]
                ]
                
                # Create new record with trimmed sequence
                new_record = SeqRecord(
                    Seq(trimmed_seq),
                    id=new_id,
                    description=f"{record.description} length={len(trimmed_seq)}",
                    letter_annotations={"phred_quality": trimmed_qual}
                )
                
                processed_records.append((new_record, 'processed'))
                
            except Exception:
                continue

        if not processed_records:
            return [(record, 'binned')]
            
        return processed_records

    def split_reads(self, input_file: str, processed_output: str, lowqual_output: str,
                   bin_output: str, stats_output: str):
        """
        Split FastQ reads based on primer pairs and indexes.
        
        Args:
            input_file: Input FastQ file path
            processed_output: Output path for processed reads
            lowqual_output: Output path for low quality reads
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
        
        processed_records = []
        lowqual_records = []
        binned_records = []
        
        with gzip.open(input_file, 'rt') as handle:
            for record in SeqIO.parse(handle, 'fastq'):
                stats['total_sequences'] += 1
                processed_segments = self.process_record(record)
                stats['total_segments'] += len(processed_segments)
                
                for new_record, classification in processed_segments:
                    if classification == 'processed':
                        processed_records.append(new_record)
                        stats['processed'] += 1
                    elif classification == 'lowqual':
                        lowqual_records.append(new_record)
                        stats['lowqual'] += 1
                    else:
                        binned_records.append(new_record)
                        stats['binned'] += 1

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

def parse_indexes(index_args: List[str]) -> Dict[str, List[str]]:
    """Parse index arguments into a dictionary."""
    index_dict = {}
    for arg in index_args:
        label, sequences = arg.split(':')
        start_seq, end_seq = sequences.split(',')
        index_dict[label] = [start_seq, end_seq]
    return index_dict

def main():
    default_indexes = {
        '1': ['AAATTTGGGCCC', 'GGGCCCAAATTT'],
        '2': ['TTTCCCAAAGGG', 'CCCTTTGGGAAA'],
        '3': ['GGGAAACCCTTT', 'AAAGGGTTTCCC'],
        '4': ['CCCGGGTTTAAA', 'TTTAAACCCGGG'],
        '5': ['AAACCCGGGAAA', 'TTTCCCGGGTTT'],
        '6': ['TTTGGGAAATTT', 'AAATTTCCCAAA'],
        '7': ['GGGTTTCCCGGG', 'CCCGGGAAACCC'],
        '8': ['CCCAAATTTCCC', 'GGGAAATTTGGG'],
        '9': ['AAAGGGAAAGGG', 'CCCTTTCCCTTT'],
        '10': ['TTTAAATTTAAA', 'TTTAAATTTAAA'],
        '11': ['GGGCCCGGGCCC', 'GGGCCCGGGCCC'],
        '12': ['CCCTTTCCCTTT', 'CCCTTTCCCTTT']
    }

    parser = argparse.ArgumentParser(description="Split FASTQ files based on index and primer sequences.")
    parser.add_argument("input_file", help="Input FASTQ file (gzipped)")
    parser.add_argument("--processed-output", required=True,
                       help="Output file for processed reads")
    parser.add_argument("--lowqual-output", required=True,
                       help="Output file for low quality reads")
    parser.add_argument("--bin-output", required=True,
                       help="Output file for binned reads")
    parser.add_argument("--stats-output", required=True,
                       help="Output statistics file")
    parser.add_argument("-fp", "--forward-primer",
                       default="AAGCAGTGGTATCAACGCAGAGT",
                       help="Forward primer sequence")
    parser.add_argument("-rp", "--reverse-primer",
                       default="ACTCTGCGTTGATACCACTGCTT",
                       help="Reverse primer sequence")
    parser.add_argument("-i", "--indexes", nargs='+',
                       help="Index sequences as 'index:start_seq,end_seq'")

    args = parser.parse_args()

    index_dict = default_indexes if args.indexes is None else parse_indexes(args.indexes)
    
    splitter = FastqSplitter(args.forward_primer, args.reverse_primer, index_dict)
    splitter.split_reads(
        args.input_file,
        args.processed_output,
        args.lowqual_output,
        args.bin_output,
        args.stats_output
    )

if __name__ == "__main__":
    main()
