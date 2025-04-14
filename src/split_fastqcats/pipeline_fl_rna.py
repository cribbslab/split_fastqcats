"""
===========================
Pipeline fl_rna
===========================

Overview
========

The aim of this pipeline is to take a nanopore input fastq and then process
the file so that it splits the files out into individual fastq files based
on the barcode sequences. It will then identify full length reads based on
appropriate read segments flanked by primers in the correct orientation and
containing polyA tails. Reads with reverse orientation are re-oriented by 
reverse complementation. Concatenated reads, identified as those with more
than two primer hits, are de-concatenated and re-oriented.
 

Configuration
-------------

The pipeline requires a configured :file:`pipeline.yml` file.
CGATReport report requires a :file:`conf.py` and optionally a
:file:`cgatreport.ini` file (see :ref:`PipelineReporting`).


To generate the default config file to change the running of the pipeline you need to
run:

split-fastqcats fl_rna config

This will generate a pipeline.yml file that the user can modify to change the
output of the pipeline.

Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning` on general
information how to use CGAT pipelines.

Once the user has modified the pipeline.yml file the
pipeline can then be ran using the following commandline command:

split-fastqcats fl_rna make full -v5

You can run the pipeline locally (without a cluster) using --local

split-fastqcats fl_rna make full -v5 --local

Input files
-----------

fastq.gz file of nanopore reads that have been sequenced with 
Macosko's TSO and Trimer RT primers. 
The pipeline.yml can be customised to allow for other primer.
Primers need to specified in the orientation expected in the sense
strand/second strand of cDNA in forward orientation. See the 
README.md for full orientation of reads expected

Pipeline output
===============

Individual fastq files split based on the presence of the primer
pairs in the correct order. Only de-concatenated reads between 
300 to 50,000 bp are kept.

Code
====

"""
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
    statement = ''' python %(PYTHON_ROOT)s/fastq_splitter_fuzzy.py -e 0.3 --num_workers 4 \
                   --processed-output %(name)s.processed.fastq.gz \
                   --lowqual-output %(name)s.lowqual.fastq.gz \
                   --bin-output %(name)s.binned.fastq.gz \
                   --stats-output %(name)s.stats.csv \
                   -res %(results_dir)s -i %(infile)s -fp %(FP)s -rp %(RP)s -v'''
    P.run(statement, job_options='-t 01:30:00', job_memory="20G", job_threads=4, without_cluster = False)

@follows(separate_by_primer_pairs)
@follows(mkdir("merged_results.dir"))
@transform('data.dir/*.fastq.gz',
           regex('data.dir/(\S+).fastq.gz'),
           r"merged_results.dir/\1.merge_complete")

def merge_by_sample(infile, outfile):
    '''
    Merge binned, lowqual, and processed fastq files from separate_samples.dir.
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
    Merge stats from all separate_samples.dir into a single CSV file.
    '''
    name = os.path.basename(infile).replace('.fastq.gz', '.stats.csv')
    PYTHON_ROOT = os.path.join(os.path.dirname(__file__), "python/")
    statement = '''python %(PYTHON_ROOT)s/merge_primer_stats.py --output merged_results.dir/%(name)s --input-dir separate_samples.dir'''
    P.run(statement)

@follows(merge_stats)
def full():
    pass

def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
