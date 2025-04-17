"""
split_fastqcats - A tool for processing and splitting FastQ reads based on primer sequences.
Helper utilities for split_fastqcats pipelines

"""

from ..version import __version__


from .fastq_splitter_by_primer import FastqSplitter as PrimerSplitter
from .fastq_splitter_by_index import FastqSplitter as IndexSplitter

__all__ = ['PrimerSplitter', 'IndexSplitter']
