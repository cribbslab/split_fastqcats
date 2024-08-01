import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import gzip

def split_reads(input_file, indexes, primer, output_file):
    with gzip.open(input_file, 'rt') as handle:
        records = list(SeqIO.parse(handle, 'fastq'))

    new_records = []
    for record in records:
        sequence = str(record.seq)
        was_split = False
        
        # Search each index and the primer in the entire read
        for index, label in indexes.items():
            search_pattern = index + primer
            pos = 0
            # Find all occurrences of the pattern
            while True:
                start_pos = sequence.find(search_pattern, pos)
                if start_pos == -1:
                    break  # No more patterns found
                pos = start_pos + len(search_pattern)  # Update position to continue search
                new_seq = sequence[start_pos + len(search_pattern):]
                new_qual = record.letter_annotations["phred_quality"][start_pos + len(search_pattern):]
                
                new_record = SeqRecord(Seq(new_seq),
                                       id=f"{record.id}_{label}",
                                       description="",
                                       letter_annotations={"phred_quality": new_qual})
                new_records.append(new_record)
                was_split = True

        # If no index-primer match was found, append the original read
        if not was_split:
            new_record = SeqRecord(Seq(sequence),
                                   id=f"{record.id}_original",
                                   description=record.description,
                                   letter_annotations={"phred_quality": record.letter_annotations["phred_quality"]})
            new_records.append(new_record)
            print(f"No match found for record ID {record.id}, appending original read.")

    with open(output_file, 'w') as output_handle:
        SeqIO.write(new_records, output_handle, 'fastq')

    print(f"Generated {len(new_records)} new reads, including originals for unmatched records.")

def parse_indexes(index_strings):
    index_dict = {}
    for item in index_strings:
        index, label = item.split(':')
        index_dict[index] = label
    return index_dict

def main():
    default_indexes = {
        'AAATTTGGGCCC': '1',
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
        'CCCTTTCCCTTT': '12'
    }

    parser = argparse.ArgumentParser(description="Split FASTQ files based on index and primer sequences.")
    parser.add_argument("input_file", type=str, help="Input FASTQ file (gzipped).")
    parser.add_argument("output_file", type=str, help="Output FASTQ file.")
    parser.add_argument("-p", "--primer", type=str, default="AAGCAGTGGTATCAACGCAGAGT", help="Primer sequence to look for. Default is 'AAGCAGTGGTATCAACGCAGAGT'.")
    parser.add_argument("-i", "--indexes", type=str, nargs='+', help="Index sequences and their labels as 'index:label'. Default is a preset list.")

    args = parser.parse_args()

    # Use default indexes if none are provided
    index_dict = default_indexes if args.indexes is None else parse_indexes(args.indexes)

    split_reads(args.input_file, index_dict, args.primer, args.output_file)

if __name__ == "__main__":
    main()
