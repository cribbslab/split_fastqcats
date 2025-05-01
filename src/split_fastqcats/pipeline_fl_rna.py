"""
===========================
Pipeline fl_rna
===========================

Full-Length RNA (fl_rna) Pipeline for split_fastqcats
-----------------------------------------------------

This pipeline is part of the split_fastqcats toolkit and is designed to process long-read (e.g., Nanopore) cDNA FastQ files, specifically for workflows where full-length RNA sequencing is performed and multiple transcript copies may be concatenated in a single read.

**Key Features:**
- Splits raw FastQ files into individual reads based on barcode and primer sequences.
- Identifies full-length reads by detecting segments flanked by expected primers in correct orientation, with polyA tail validation.
- Automatically re-orients reads with reverse orientation.
- De-concatenates reads containing multiple transcript copies (multiple primer pairs in a single read).
- Merges and summarizes results across samples, providing detailed QC statistics.
- Supports high-throughput batch processing for both local and cluster environments.

**Inputs:**
- Nanopore (or similar) FastQ files with concatenated cDNA reads
- User-customizable pipeline.yml configuration file

**Outputs:**
- Processed, low-quality, and binned FastQ files for each sample
- Merged statistics and QC summaries

**Usage:**
1. Generate a config file:
   ```bash
   split-fastqcats fl_rna config
   ```
2. Edit `pipeline.yml` as needed
3. Run the pipeline:
   ```bash
   split-fastqcats fl_rna make full -v5
   # or locally
   split-fastqcats fl_rna make full -v5 --local
   ```

This pipeline leverages CGAT-style batch processing and is tightly integrated with the split_fastqcats core logic for robust, error-tolerant splitting and de-concatenation of long-read RNA data.
"""

# ===================
# Input files
# ===================
#
# fastq.gz file of nanopore reads that have been sequenced with 
# Macosko's TSO and Trimer RT primers. 
# The pipeline.yml can be customised to allow for other primer.
# Primers need to specified in the orientation expected in the sense
# strand/second strand of cDNA in forward orientation. See the 
# README.md for full orientation of reads expected
#
# ===================
# Pipeline output
# ===================
#
# Individual fastq files split based on the presence of the primer
# pairs in the correct order. Only de-concatenated reads between 
# 300 to 50,000 bp are kept.
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

# Get the list of indexes from the YAML config file

SEQUENCESUFFIXES = ("*.fastq.gz")

FASTQTARGET = tuple([os.path.join("data.dir/", suffix_name)
                       for suffix_name in SEQUENCESUFFIXES])

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
    statement = '''zcat %(infile)s | split -l %(split)s --additional-suffix=.fastq - %(name)s. && mv %(name)s*.fastq split_tmp.dir/ '''
    P.run(statement)

@follows(split_fastq)
@follows(mkdir("separate_samples.dir"))
@transform('split_tmp.dir/*.fastq',
           regex("split_tmp.dir/(\S+).fastq"),
           r"separate_samples.dir/\1/\1.stats.csv")

def separate_by_primer_pairs(infile, outfile):
    '''
    Identify barcode and split reads accordingly using a different script.
    '''
    PYTHON_ROOT = os.path.join(os.path.dirname(__file__), "python/")
    name = os.path.basename(infile).replace('.fastq', '')
    results_dir = os.path.join("separate_samples.dir", name)
    statement = ''' python %(PYTHON_ROOT)s/fastq_splitter_by_primer.py -e %(error)s --num_workers 4 \
                   --processed-output %(name)s.processed.fastq.gz \
                   --lowqual-output %(name)s.lowqual.fastq.gz \
                   --bin-output %(name)s.binned.fastq.gz \
                   --stats-output %(name)s.stats.csv \
                   -res %(results_dir)s -i %(infile)s -fp %(FP)s -rp %(RP)s -v'''
    P.run(statement, job_memory=PARAMS.get('splitter_mem'), job_threads=4, without_cluster = False)

@follows(separate_by_primer_pairs)
@follows(mkdir("merged_results.dir"))
@transform('data.dir/*.fastq.gz',
           regex('data.dir/(\S+).fastq.gz'),
           r"merged_results.dir/\1.merge_complete")

def merge_by_sample(infile, outfile):
    '''
    Merge all binned, low-quality, and processed FastQ files for a given sample from separate_samples.dir
    into consolidated files in merged_results.dir. This step ensures all read categories are available
    for downstream analysis and reporting.
    '''
    name = os.path.basename(infile).replace('.fastq.gz', '')
    statement = '''cat separate_samples.dir/%(name)s.*/*.binned.fastq.gz > merged_results.dir/%(name)s.binned.fastq.gz &&
                   cat separate_samples.dir/%(name)s.*/*.lowqual.fastq.gz > merged_results.dir/%(name)s.lowqual.fastq.gz &&
                   cat separate_samples.dir/%(name)s.*/*.processed.fastq.gz > merged_results.dir/%(name)s.processed.fastq.gz &&
                   touch %(outfile)s'''
    P.run(statement)

@follows(merge_by_sample)
@transform('data.dir/*.fastq.gz',
           regex('data.dir/(\S+).fastq.gz'),
           r"merged_results.dir/\1.merged_stats.csv")
def merge_stats(infile, outfile):
    '''
    Merge all statistics CSV files from separate_samples.dir into a single summary CSV file in merged_results.dir.
    This provides a unified QC and processing summary for each sample.
    '''
    name = os.path.basename(infile).replace('.fastq.gz', '.stats.csv')
    PYTHON_ROOT = os.path.join(os.path.dirname(__file__), "python/")
    statement = '''python %(PYTHON_ROOT)s/merge_primer_stats.py --output merged_results.dir/%(name)s --input-dir separate_samples.dir'''
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
