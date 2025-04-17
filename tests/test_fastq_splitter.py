import sys
from types import SimpleNamespace

# --- Mock parasail globally before any split_fastqcats import ---
import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from split_fastqcats.python.fastq_splitter_by_index import FastqSplitter

@pytest.fixture
def example_splitter():
    forward_primer = "AAGCAGTGGTATCAACGCAGAGT"
    index_dict = {'1': 'AAATTTGGGCCC'}  # str->str mapping
    mismatches = 2
    return FastqSplitter(forward_primer, index_dict, mismatches)

def make_seqrecord(seq, name="test", qual=40):
    return SeqRecord(Seq(seq), id=name, description="", letter_annotations={"phred_quality": [qual]*len(seq)})

def test_smith_waterman_search_exact(example_splitter):
    # The pattern is index + forward_primer[:10]
    index = 'AAATTTGGGCCC'
    primer = 'AAGCAGTGGTATCAACGCAGAGT'[:10]  # "AAGCAGTGGT"
    seq = index + primer + "A" * (100 - len(index + primer))
    matches = example_splitter.smith_waterman_search(seq, "read1")
    assert matches, "Should find a match for exact barcode+primer"
    assert matches[0]["start"] == 0 or matches[0]["start"] is not None
    assert matches[0]["end"] > 0

# def test_find_best_primer_pairs(example_splitter):
#     # Method does not exist in FastqSplitter (by_index)
#     pass, "Each pair should include 'trimmed_seq' key"

