# split_fastqcats

A tool for processing and splitting concatenated FastQ reads based on primer sequences, particularly designed for long-read RNA sequencing data.

## Overview

`split_fastqcats` is designed to handle FastQ reads with the following structure:

```
TSO---cDNA---polyA---UMI---revRTprimer
```

Example structure:
```
AAGCAGTGGTATCAACGCAGAGTGAAT---cDNA---polyA--N---UMI---GTACTCTGCGTTGATACCACTGCTT
```

The tool performs the following:
- Identifies and validates primer sequences using Smith-Waterman alignment
- Handles potential sequencing errors with configurable mismatch tolerance
- Processes polyA tails and UMIs
- Splits concatenated reads into individual segments
- Provides detailed statistics about the processing

## Installation

### Using pip

```bash
pip install split_fastqcats
```

### From source

```bash
git clone https://github.com/cribbslab/split_fastqcats.git
cd split_fastqcats
pip install -e .
```

### Using conda/mamba

```bash
mamba env create -f conda/environments/split_fastqcats.yml
conda activate fastcats
pip install -e .
```

## Usage

### Command Line Interface

Basic usage:
```bash
split-fastqcats input.fastq.gz \
    -fp AAGCAGTGGTATCAACGCAGAGT \
    -rp ACTCTGCGTTGATACCACTGCTT \
    processed.fastq.gz \
    lowqual.fastq.gz \
    bin.fastq.gz \
    stats.csv
```

Arguments:
- `input.fastq.gz`: Input FastQ file (gzipped)
- `-fp, --forward_primer`: Forward primer sequence (default: AAGCAGTGGTATCAACGCAGAGT)
- `-rp, --reverse_primer`: Reverse primer sequence (default: ACTCTGCGTTGATACCACTGCTT)
- `-i, --indexes`: Optional index sequences as 'index:start_seq,end_seq'
- `processed.fastq.gz`: Output file for correctly processed reads
- `lowqual.fastq.gz`: Output file for low quality reads
- `bin.fastq.gz`: Output file for binned reads
- `stats.csv`: Output statistics file

### Python API

```python
from split_fastqcats import FastqSplitter

# Initialize the splitter
splitter = FastqSplitter(
    forward_primer="AAGCAGTGGTATCAACGCAGAGT",
    reverse_primer="ACTCTGCGTTGATACCACTGCTT",
    index_dict={'1': ['AAATTTGGGCCC', 'GGGCCCAAATTT']}
)

# Process reads
splitter.split_reads(
    input_file="input.fastq.gz",
    processed_output="processed.fastq.gz",
    lowqual_output="lowqual.fastq.gz",
    bin_output="bin.fastq.gz",
    stats_output="stats.csv"
)
```

## Output Files

1. `processed.fastq.gz`: Contains reads that:
   - Have valid primer pairs
   - Meet quality thresholds
   - Have correct index sequences (if specified)

2. `lowqual.fastq.gz`: Contains reads that:
   - Have only one valid primer
   - Have mismatched indexes
   - Meet basic quality requirements but fail stricter criteria

3. `bin.fastq.gz`: Contains reads that:
   - Lack valid primers
   - Have too many primer matches (>10)
   - Fail basic quality requirements

4. `stats.csv`: Contains processing statistics:
   - Total sequences processed
   - Number of processed reads
   - Number of low-quality reads
   - Number of binned reads

## Quality Control Parameters

The tool implements several QC measures:
- Smith-Waterman alignment for primer detection
- Maximum 5 mismatches allowed in primer sequences
- Minimum 75% match score required for primer identification
- Index validation with 83% minimum match score
- Position validation for index sequences

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Citation

If you use this tool in your research, please cite:

```bibtex
@software{split_fastqcats2024,
  author = {Cribbs, Adam},
  title = {split_fastqcats: A tool for processing concatenated FastQ reads},
  year = {2024},
  publisher = {GitHub},
  url = {https://github.com/cribbslab/split_fastqcats}
}
