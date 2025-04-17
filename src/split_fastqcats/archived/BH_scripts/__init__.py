"""
split_fastqcats - A tool for processing and splitting FastQ reads based on primer sequences
"""

__version__ = '0.1.0'

from .fastq_splitter import FastqSplitter, main

__all__ = ['FastqSplitter', 'main']