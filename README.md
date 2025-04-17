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
- Processes polyA tails and polyT tails in reverse oriented reads
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

# Primer Pair Split
split-fastqcats primer_pair_split \
    -i input.fastq.gz \
    --processed-output processed.fastq.gz \
    --lowqual-output lowqual.fastq.gz \
    --bin-output bin.fastq.gz \
    --stats-output stats.csv \
    -fp AAGCAGTGGTATCAACGCAGAGTGAAT \
    -rp GTACTCTGCGTTGATACCACTGCTT \ ## Note primer input orientation
    --error 0.3 \
    --chunk-size 1000 \
    --num_workers 4 \
    -v

# Barcode Split
split-fastqcats barcode_split \
    -i input.fastq.gz \
    --indexes AAATTTGGGCCC TTTCCCAAAGGG \
    --processed-output processed.fastq.gz \
    --lowqual-output lowqual.fastq.gz \
    --bin-output bin.fastq.gz \
    --stats-output stats.csv \
    -fp AAGCAGTGGTATCAACGCAGAGT \
    --error 3 \
    --chunk-size 1000 \
    --num_workers 4 \
    -v

```
## Arguments

### Primer Pair Split

- `-i, --input_file`: Input FASTQ file (gzipped)
- `-res, --results-dir`: Output directory (default is based on input filename)
- `--processed-output`: Output file for correctly processed reads
- `--lowqual-output`: Output file for low-quality reads
- `--bin-output`: Output file for binned reads
- `--stats-output`: Output statistics file (CSV)
- `-fp, --forward-primer`: Forward primer sequence (default: `AAGCAGTGGTATCAACGCAGAGTGAAT`)
- `-rp, --reverse-primer`: Reverse primer sequence (default: `GTACTCTGCGTTGATACCACTGCTT`)
- `-e, --error`: Number of allowed mismatches (default: `0.3`)
- `--chunk-size`: Number of reads per chunk (default: `1000`)
- `--num_workers`: Number of parallel workers (default: `4`)
- `-v, --verbose`: Enable detailed logging

### Barcode/Index Split

- `-i, --input_file`: Input FASTQ file (gzipped)
- `-res, --results-dir`: Output directory (default is based on input filename)
- `--processed-output`: Output file for processed reads, one per index
- `--lowqual-output`: Output file for low-quality reads
- `--bin-output`: Output file for binned reads (no barcode)
- `--stats-output`: Output statistics file (CSV)
- `-fp, --forward-primer`: Forward primer sequence (default: `AAGCAGTGGTATCAACGCAGAGT`)
- `--indexes`: List of index/barcode sequences
- `-e, --error`: Number of allowed mismatches (default: `3`)
- `--chunk-size`: Number of reads per chunk (default: `1000`)
- `--num_workers`: Number of parallel workers (default: `4`)
- `-v, --verbose`: Enable detailed logging

---

## Outputs

### Primer Pair Split Outputs

1. **`processed.fastq.gz`**: Contains reads that:
   - Have valid primer pairs
   - Meet quality thresholds
   - Have correct index sequences (if specified)

2. **`lowqual.fastq.gz`**: Contains reads that:
   - Have only one valid primer
   - Have mismatched indexes
   - Meet basic quality requirements but fail stricter criteria (e.g., length, polyA)

3. **`bin.fastq.gz`**: Contains reads that:
   - Lack valid primers
   - Have too many primer matches (>10)
   - Fail basic quality requirements

4. **`stats.csv`**: Contains processing statistics:
   - Total sequences processed
   - Number of processed reads
   - Number of low-quality reads
   - Number of binned reads
   - Full-length vs low-quality segment breakdown

---

### Barcode/Index Split Outputs

1. **`processed_index_<barcode>.fastq.gz`**: One output file per index containing processed reads for that barcode.

2. **`lowqual.fastq.gz`**: Contains low-quality reads that fail length filters.

3. **`bin.fastq.gz`**: Contains reads with no barcode hits.

4. **`stats.csv`**: Contains processing statistics:
   - Total sequences processed
   - Number of processed reads
   - Number of binned reads (no barcode)
   - Total segments identified by de-concatenation
   - Breakdown of processed vs low-quality segments failing on length filter
   - Barcode Count Table: A separate table giving a count of segments for each barcode in the list.


The tool supports CGAT-style pipelines for batch processing. 
Raw fastq.gz files to be placed is directory in the working directory called data.dir. 

```bash
# To de-concatenate and de-multiplex reads by index/barcode on a cluster (e.g. Slurmm). 
split-fastqcats split_by_index config         # Generate pipeline.yml
split-fastqcats split_by_index make full -v5  # Run the pipeline


## To identify full length reads on a cluster (e.g. Slurmm)
split-fastqcats fl_rna config         # Generate pipeline.yml
split-fastqcats fl_rna make full -v5  # Run the pipeline


## To run the pipelines locally without a cluster:
split-fastqcats split_by_index make full -v5 --local
split-fastqcats fl_rna make full -v5 --local

```
## Outputs for each sample will be produced in the merged_results.dir.

### Python API

```python

# To de-concatenate and identify full length reads
from split_fastqcats import PrimerSplitter

splitter = PrimerSplitter(
    forward_primer="AAGCAGTGGTATCAACGCAGAGTGAAT",
    reverse_primer="GTACTCTGCGTTGATACCACTGCTT",
    error=0.3
)

splitter.parallel_split_reads(
    input_file="input.fastq.gz",
    processed_output="processed.fastq.gz",
    lowqual_output="lowqual.fastq.gz",
    bin_output="bin.fastq.gz",
    stats_output="stats.csv",
    num_workers=4,
    chunk_size=1000
)

# To de-concatenate and de-multiplex by barcodes
from split_fastqcats import IndexSplitter

index_dict = {
    '1': 'AAATTTGGGCCC',
    '2': 'TTTCCCAAAGGG'
}

splitter = IndexSplitter(
    forward_primer="AAGCAGTGGTATCAACGCAGAGT",
    index_dict=index_dict,
    error=3
)

splitter.parallel_split_reads(
    input_file="input.fastq.gz",
    processed_output="processed.fastq.gz",
    lowqual_output="lowqual.fastq.gz",
    bin_output="bin.fastq.gz",
    stats_output="stats.csv",
    num_workers=4,
    chunk_size=1000
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

The primer pair matching de-concatenating tool used to identify full-length mRNA implements several QC measures:
- Smith-Waterman alignment for primer detection
- Minimum 70% match score required for primer identification (using default parameters)
- Minimum sequence length 300, Max length 50,000 bp
- Checking for polyA tails at the ends of the sequences

The barcode matching de-concatenating tool uses:
- Index validation with 80% minimum match score (using default parameters)
- Position validation for index sequences - i.e. they should flank the primers

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
