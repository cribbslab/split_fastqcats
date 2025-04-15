"""
split_fastqcats - A tool for processing and splitting FastQ reads based on primer sequences.
Helper utilities for split_fastqcats pipelines

"""

__version__ = '0.0.1'

from .fastq_splitter_fuzzy import FastqSplitter as PrimerSplitter
from .fastq_splitter_index_fuzzy import FastqSplitter as IndexSplitter

__all__ = ['PrimerSplitter', 'IndexSplitter']