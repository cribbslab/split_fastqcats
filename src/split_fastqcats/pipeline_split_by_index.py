"""
===========================
Pipeline split_by_index
===========================

Barcode/Index-Based Splitting Pipeline for split_fastqcats
----------------------------------------------------------

This pipeline is part of the split_fastqcats toolkit and is designed to process long-read (e.g., Nanopore) cDNA FastQ files where demultiplexing and de-concatenation are required based on barcode (index) sequences.

**Key Features:**
- Demultiplexes and splits concatenated reads into individual FastQ files based on user-defined barcode/index sequences.
- Uses robust, error-tolerant matching to handle sequencing errors and protocol variability.
- Supports batch processing for both local and cluster (e.g., Slurm) environments.
- Merges and summarizes results across samples, providing detailed QC statistics.

**Inputs:**
- Nanopore (or similar) FastQ files with barcoded, concatenated cDNA reads
- User-customizable pipeline.yml configuration file specifying index/barcode sequences

**Outputs:**
- Processed, low-quality, and binned FastQ files for each barcode/sample
- Merged statistics and QC summaries

**Usage:**
1. Generate a config file:
   ```bash
   split-fastqcats split_by_index config
   ```
2. Edit `pipeline.yml` to specify barcode/index sequences and other options
3. Run the pipeline:
   ```bash
   split-fastqcats split_by_index make full -v5
   # or locally
   split-fastqcats split_by_index make full -v5 --local
   ```

This pipeline leverages CGAT-style batch processing and is tightly integrated with the split_fastqcats core logic for robust, error-tolerant demultiplexing and de-concatenation of long-read RNA data.
"""

# ===================
# Pipeline output
# ===================
#
# Individual fastq files split based on the presence of the barcode.
# Only de-concatenated reads between 50 to 50,000 bp are kept. 
#
# ===================
# Code
# ===================

import sys
import os
import glob
import pandas as pd
from ruffus import *
import cgatcore.pipeline as P
import cgatcore.experiment as E

# Load parameters from the config file
PARAMS = P.get_parameters([
    "%s/pipeline.yml" % os.path.splitext(__file__)[0],
    "../pipeline.yml",
    "pipeline.yml"])


SEQUENCESUFFIXES = ("*.fastq.gz")

FASTQTARGET = tuple([os.path.join("data.dir/", suffix_name)
                       for suffix_name in SEQUENCESUFFIXES])

# Get the list of indexes from the YAML config file
INDEXES = PARAMS.get('indexes', [])
PYTHON_ROOT = os.path.join(os.path.dirname(__file__), "python/")
BASH_ROOT = os.path.join(os.path.dirname(__file__), "bash/")

@follows(mkdir("split_tmp.dir"))
@transform('data.dir/*.fastq.gz',
           regex('data.dir/(\S+).fastq.gz'),
           r"split_tmp.dir/\1.aa.fastq")
def split_fastq(infile, outfile):
    '''
    Split the fastq file into smaller chunks.
    '''
    infile = "".join(infile)
    name = infile.replace('data.dir/','').replace('.fastq.gz','')
    statement = '''zcat %(infile)s | split -l %(split)s --additional-suffix=.fastq - %(name)s. &&
                   mv %(name)s*.fastq split_tmp.dir/'''
    P.run(statement)

@follows(split_fastq)
@follows(mkdir("separate_samples.dir"))
@transform('split_tmp.dir/*.fastq',
           regex("split_tmp.dir/(\S+).fastq"),
           r"separate_samples.dir/\1/\1.stats.csv")
def separate_by_index(infile, outfile):
    '''
    Identify barcode and split reads accordingly using a different script.
    '''
    PYTHON_ROOT = os.path.join(os.path.dirname(__file__), "python/")
    name = os.path.basename(infile).replace('.fastq', '')
    results_dir = os.path.join("separate_samples.dir", name)
    barcodes = " ".join(INDEXES)
    statement = '''mkdir -p %(results_dir)s &&
                   python %(PYTHON_ROOT)s/fastq_splitter_index_fuzzy.py -e %(error)s --num_workers 4 \
                   --processed-output %(name)s.processed \
                   --lowqual-output %(name)s.lowqual.fastq.gz \
                   --bin-output %(name)s.binned.fastq.gz \
                   --stats-output %(name)s.stats.csv \
                   -res %(results_dir)s -i %(infile)s -fp %(FP)s --indexes %(barcodes)s -v'''
    P.run(statement, job_options='-t 03:00:00', job_memory="20G", job_threads=4, without_cluster = False)

@follows(separate_by_index)
@follows(mkdir("merged_results.dir"))
@transform('data.dir/*.fastq.gz',
           regex('data.dir/(\S+).fastq.gz'),
           r"merged_results.dir/\1.merge_complete")
def merge_by_index(infile, outfile):
    '''
    Merge all binned, low-quality, and processed FastQ files for a given barcode/sample from separate_samples.dir
    into consolidated files in merged_results.dir. This step ensures all read categories are available
    for downstream analysis and reporting.
    '''
    BASH_ROOT = os.path.join(os.path.dirname(__file__), "bash/")
    barcodes = " ".join(INDEXES)
    name = os.path.basename(infile).replace('.fastq.gz', '')
    statement = '''bash %(BASH_ROOT)s/cat_files_by_index.sh "%(barcodes)s" "%(name)s" && touch %(outfile)s'''
    P.run(statement)

@follows(merge_by_index)
@transform('data.dir/*.fastq.gz',
           regex('data.dir/(\S+).fastq.gz'),
           r"merged_results.dir/\1.index_stats.csv")
def merge_stats(infile, outfile):
    '''
    Merge all statistics CSV files from separate_samples.dir into a single summary CSV file in merged_results.dir.
    This provides a unified QC and processing summary for each barcode/sample.
    '''
    PYTHON_ROOT = os.path.join(os.path.dirname(__file__), "python/")
    name = os.path.basename(infile).replace('.fastq.gz', '')
    statement = '''python %(PYTHON_ROOT)s/merge_index_stats.py --input-dir separate_samples.dir --output-file merged_results.dir/%(name)s.index_stats.csv'''
    P.run(statement)

@follows(merge_stats)
def full():
    '''
    Dummy target for the pipeline. Ensures all steps (splitting, separating, merging, summarizing)
    are executed in the correct order when running the full workflow.
    '''
    pass

def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
