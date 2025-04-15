##Edits to Adam's main 'fastq_splitter.py' script 
# - Reverse complement polyT sequences 
# - Reintroduced UMIs for QC processing
# - Filters reads without polyA
# - Parallelised for efficient processing of large patient samples

import regex
import sys
import argparse
import gzip
from Bio import SeqIO, Align
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import csv
from typing import Dict, List, Tuple, Optional
from multiprocessing import Pool, cpu_count

class FastqSplitter:
    def __init__(self, forward_primer: str, reverse_primer: str, index_dict: Dict[List[str], str]):
        """
        Initialize the FastQ splitter with primers and index dictionary.
        
        Args:
            forward_primer: Forward primer sequence
            reverse_primer: Reverse primer sequence
            index_dict: Dictionary mapping index numbers 
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
                    trimmed_seq = sequence[pos1_end:pos2_start]
                elif type1 == 'reverse' and type2 == 'forward':
                    # Reverse complement forward->reverse order - adapted from Adam's script to not just swap the foward and reverse but the whole seq
                    match1, match2 = match2, match1
                    pos1_start, pos2_start = pos2_start, pos1_start
                    pos1_end, pos2_end = pos2_end, pos1_end
                    valid_pair = True
                    trimmed_seq = str(Seq(sequence[pos1_end:pos2_start]).reverse_complement())
                
                if valid_pair:
                    # Check if the distance between primers is reasonable
                    distance = pos2_start - pos1_end
                    if 10 <= distance <= 2000:
                        segment = {
                            'trimmed_seq': trimmed_seq,
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

    def process_record(self, records: List[SeqRecord]) -> List[Tuple[SeqRecord, str]]:
        """
        Process a single FastQ record and split it into multiple records if needed.
        
        Args:
            record: SeqRecord object to process
            
        Returns:
            List of tuples (processed record, classification)
        """
        for record in records:
            seq = str(record.seq)
            matches = self.find_best_primer_pairs(seq)
        
            if not matches:
                return [(record, 'binned')]

            processed_records = []
            for i, match in enumerate(matches, 1):
                check_trimmed_seq = match['trimmed_seq']
            
                # Skip segments that are too short or too long
                if len(check_trimmed_seq) < 50 or len(check_trimmed_seq) > 5000:
                    continue
                
                umi_candidate = check_trimmed_seq[-12:]
                umi_pattern = regex.compile(r'((?:(?:A{3}|T{3}|C{3}|G{3})){4}){i<=0,e<=1}')
                umi_match = umi_pattern.fullmatch(umi_candidate)
            
                umi_seq = umi_match.group(1) if umi_match else None
                new_id = f"{record.id}_segment_{i}" + (f"_{umi_seq}" if umi_seq else "")
            
                polyA = regex.findall("(AAAAAAAA){e<=0}", str(seq))
            
                if umi_seq and polyA:  
                    classification = "processed"  # UMI and polyA present = processed
                else:
                    classification = "lowqual"   # Either missing = lowqual
                
                try:
                    #Get quality of trimmed seq
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
                
                    results_records.append((new_record, classification))
                
                except Exception:
                    continue

        if not results:
            return [(record, 'binned')]
            
        return results_records

    
    def worker(self, args: str):
        records_chunk = args
        return self.process_record(records_chunk)

    def parallel_split_reads(self, input_file: str, processed_output: str, lowqual_output: str,
                   bin_output: str, stats_output: str):
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
        
        processed_records = []
        lowqual_records = []
        binned_records = []
        
        with gzip.open(input_file, 'rt') as handle:
            records = list(SeqIO.parse(handle, 'fastq'))
            for record in records:
                stats['total_sequences'] += 1
            # Split records into chunks
            num_workers = cpu_count()
            chunk_size = len(records) // num_workers
            chunks = [records[i:i + chunk_size] for i in range(0, len(records), chunk_size)]
    
            # Prepare arguments for worker processes
            args = [(chunk) for chunk in chunks]

            # Parallel processing
            with Pool(num_workers) as pool:
                results = pool.map(self.worker, args)

            all_results = [record_tuple for sublist in results for record_tuple in sublist]
	    stats['total_segments'] = len(all_results)
                
            for new_record, classification in all_results:
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

def parse_indexes(index_args: List[str]) -> Dict[List[str], str]:
    """Parse index arguments into a dictionary."""
    index_dict = {}
    for arg in index_args:
        sequences, label = arg.split(':')
    return index_dict

def main():
    #Indexes corrected: should not be reverse transcribed? eg. no need for complement as only one seq
    default_indexes = {'AAATTTGGGCCC': '1',
                       'TTTCCCAAAGGG': '2',
                       'GGGAAACCCTTT': '3',
                       'CCCGGGTTTAAA': '4',
                       'AAACCCGGGAAA': '5',
                       'TTTGGGAAATTT': '6',
                       'GGGTTTCCCGGG': '7',
                       'CCCAAATTTCCC': '8',
                       'AAAGGGAAAGGG': '9',
                       'TTTAAATTTAAA': '10',
                       'GGGCCCGGGCCC': '11',
                       'CCCTTTCCCTTT': '12'}

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
                       default="AAGCAGTGGTATCAACGCAGAGTGAATGGG",
                       help="Forward primer sequence")
    parser.add_argument("-rp", "--reverse-primer",
                       default="GTACTCTGCGTTGATACCACTGCTT",
                       help="Reverse primer sequence")
    parser.add_argument("-i", "--indexes", nargs='+',
                       help="Index sequences as 'index:start_seq,end_seq'")

    args = parser.parse_args()

    index_dict = default_indexes if args.indexes is None else parse_indexes(args.indexes)
    
    splitter = FastqSplitter(args.forward_primer, args.reverse_primer, index_dict)
    splitter.parallel_split_reads(
        args.input_file,
        args.processed_output,
        args.lowqual_output,
        args.bin_output,
        args.stats_output
    )

if __name__ == "__main__":
    main()
