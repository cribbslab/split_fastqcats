"""
split_fastqcats - A tool for processing and splitting FastQ reads based on primer sequences.

"""

# split_fastqcats/__init__.py
from .version import __version__

from .python import PrimerSplitter, IndexSplitter

__all__ = ['PrimerSplitter', 'IndexSplitter']
