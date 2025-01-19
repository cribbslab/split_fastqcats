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
        min_score_threshold = 0.75 * primer_length * self.aligner.match_score

        for i in range(len(sequence) - primer_length + 1):
            sub_seq = sequence[i:i + primer_length]
            alignments = self.aligner.align(sub_seq, primer)
            
            if not alignments:
                continue
                
            alignment = alignments[0]
            score = alignment.score
            
            if score >= min_score_threshold:
                mismatches = sum(1 for a, b in zip(sub_seq, primer) if a != b)
                if mismatches <= max_mismatches:
                    matches.append({
                        'start': i,
                        'end': i + primer_length,
                        'score': score,
                        'mismatches': mismatches
                    })
                    
        return sorted(matches, key=lambda x: (x['mismatches'], x['start']))

    def find_best_primer_pairs(self, sequence: str) -> List[dict]:
        """
        Find the best matching primer pairs in the sequence.
        
        Args:
            sequence: Input DNA sequence
            
        Returns:
            List of dictionaries containing paired sequence information
        """
        forward_matches = self.smith_waterman_search(sequence, self.forward_primer)
        reverse_matches = self.smith_waterman_search(sequence, self.reverse_primer)
        
        paired_sequences = []
        used_reverse_indices = set()

        for f_match in forward_matches:
            best_pair = None
            best_distance = float('inf')
            best_mismatches = float('inf')

            for r_match in reverse_matches:
                if r_match['start'] > f_match['end'] and r_match['start'] not in used_reverse_indices:
                    distance = r_match['start'] - f_match['end']
                    total_mismatches = f_match['mismatches'] + r_match['mismatches']
                    
                    if distance < best_distance or (distance == best_distance and total_mismatches < best_mismatches):
                        best_distance = distance
                        best_mismatches = total_mismatches
                        best_pair = r_match
            
            if best_pair:
                paired_sequences.append({
                    'trimmed_seq': sequence[f_match['end']:best_pair['start']],
                    'forward': (f_match['start'], f_match['end']),
                    'reverse': (best_pair['start'], best_pair['end']),
                    'total_mismatches': f_match['mismatches'] + best_pair['mismatches']
                })
                used_reverse_indices.add(best_pair['start'])

        return paired_sequences

    def process_record(self, record: SeqRecord) -> Tuple[List[SeqRecord], str]:
        """
        Process a single FastQ record.
        
        Args:
            record: SeqRecord object to process
            
        Returns:
            Tuple of (processed records, classification)
        """
        seq = str(record.seq)
        matches = self.find_best_primer_pairs(seq)
        
        if not matches or len(matches) > 10:
            return [], 'binned'

        processed_records = []
        for match in matches:
            trimmed_seq = match['trimmed_seq']
            
            # Check for indexes
            forward_found = False
            reverse_found = False
            
            for index_pair in self.index_dict.values():
                start_index, end_index = index_pair
                
                f_matches = self.smith_waterman_search(trimmed_seq, start_index)
                if f_matches and f_matches[0]['score'] >= 0.83 * len(start_index) * 2:
                    forward_found = True
                    
                r_matches = self.smith_waterman_search(trimmed_seq, end_index)
                if r_matches and r_matches[0]['score'] >= 0.83 * len(end_index) * 2:
                    reverse_found = True
                    
                if forward_found and reverse_found:
                    break

            # Create new record with trimmed sequence
            new_record = SeqRecord(
                Seq(trimmed_seq),
                id=record.id,
                description=record.description,
                letter_annotations={
                    "phred_quality": record.letter_annotations["phred_quality"][
                        match['forward'][1]:match['reverse'][0]
                    ]
                }
            )
            
            if forward_found and reverse_found:
                return [new_record], 'processed'
            elif forward_found or reverse_found:
                return [new_record], 'lowqual'
                
        return [], 'binned'

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
                records, classification = self.process_record(record)
                
                if classification == 'processed':
                    processed_records.extend(records)
                    stats['processed'] += 1
                elif classification == 'lowqual':
                    lowqual_records.extend(records)
                    stats['lowqual'] += 1
                else:
                    binned_records.append(record)
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
    parser.add_argument("-fp", "--forward_primer",
                       default="AAGCAGTGGTATCAACGCAGAGT",
                       help="Forward primer sequence")
    parser.add_argument("-rp", "--reverse_primer",
                       default="ACTCTGCGTTGATACCACTGCTT",
                       help="Reverse primer sequence")
    parser.add_argument("-i", "--indexes", nargs='+',
                       help="Index sequences as 'index:start_seq,end_seq'")
    parser.add_argument("processed_output", help="Output for processed reads")
    parser.add_argument("lowqual_output", help="Output for low quality reads")
    parser.add_argument("bin_output", help="Output for binned reads")
    parser.add_argument("stats_output", help="Output statistics file")

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
