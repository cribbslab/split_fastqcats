import sys
from types import SimpleNamespace

# --- Mock parasail globally before any split_fastqcats import ---
import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from split_fastqcats.python.fastq_splitter_by_index import FastqSplitter

@pytest.fixture
def example_splitter():
    forward_primer = "AAGCAGTGGT"
    index_dict = {"1": "AAATTTGGGCCC"}
    mismatches = 2
    return FastqSplitter(forward_primer, index_dict, mismatches)

def make_seqrecord(seq, name="test", qual=40):
    return SeqRecord(Seq(seq), id=name, description="", letter_annotations={"phred_quality": [qual]*len(seq)})

def test_smith_waterman_search_exact(example_splitter):
    index = 'AAATTTGGGCCC'
    primer = 'AAGCAGTGGTATCAACGCAGAGT'[:10]
    seq = index + primer
    matches = example_splitter.smith_waterman_search(seq, "read1")
    assert matches, "Should find a match for exact barcode+primer"
    assert matches[0]["start"] == 0 or matches[0]["start"] is not None
    assert matches[0]["end"] > 0

def test_smith_waterman_search_with_mismatch(example_splitter):
    index = 'AAATTTGGGCCC'
    primer = 'AAGCAGTGGTATCAACGCAGAGT'[:10]
    # Introduce a single mismatch in the index region
    seq = 'AAATTTGGGACC' + primer  # C->A mutation at position 10
    matches = example_splitter.smith_waterman_search(seq, "read2")
    assert matches, "Should tolerate one mismatch"
    assert matches[0]["start"] == 0 or matches[0]["start"] is not None
    assert matches[0]["end"] > 0

def test_smith_waterman_search_no_match(example_splitter):
    seq = "GGGGGGGGGGGGGGGGGGGGGGGG"
    matches = example_splitter.smith_waterman_search(seq, "read3")
    assert not matches or all(m["score"] == 0 for m in matches), "Should not find a match with wrong sequence"

def test_multiple_matches(example_splitter):
    index = 'AAATTTGGGCCC'
    primer = 'AAGCAGTGGTATCAACGCAGAGT'[:10]
    pattern = index + primer
    seq = pattern + 'NNNNN' + pattern
    matches = example_splitter.smith_waterman_search(seq, "read4")
    assert len(matches) >= 2, "Should find two matches"

def test_empty_sequence(example_splitter):
    matches = example_splitter.smith_waterman_search("", "empty")
    assert matches == [] or all(m["score"] == 0 for m in matches), "Should return empty list for empty sequence"
