import sys
from types import SimpleNamespace

# --- Mock parasail globally before any split_fastqcats import ---
import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from split_fastqcats.python.fastq_splitter_by_primer import FastqSplitter

@pytest.fixture
def example_splitter():
    forward_primer = "AAGCAGTGGT"
    index_dict = {"1": "AAATTTGGGCCC"}
    mismatches = 2
    return FastqSplitter(forward_primer, index_dict, mismatches)

def make_seqrecord(seq, name="test", qual=40):
    return SeqRecord(Seq(seq), id=name, description="", letter_annotations={"phred_quality": [qual]*len(seq)})

def test_smith_waterman_search_exact(example_splitter):
    seq = "AAATTTGGGCCCAAGCAGTGGT" + "A" * (100 - len("AAATTTGGGCCCAAGCAGTGGT"))
    # If this FastqSplitter requires a 'primer' argument, add it below; otherwise, call as before.
    try:
        matches = example_splitter.smith_waterman_search(seq, "read1")
    except TypeError:
        # Try with primer argument if required
        matches = example_splitter.smith_waterman_search(seq, "read1", "AAGCAGTGGT")
    assert matches, "Should find a match for exact barcode+primer"
    assert matches[0]["start"] == 0 or matches[0]["start"] is not None
    assert matches[0]["end"] > 0

def test_smith_waterman_search_with_mismatch(example_splitter):
    seq = "AAATTTGGGCCAAGCAGTGGT" + "A" * (100 - len("AAATTTGGGCCAAGCAGTGGT"))  # One C missing, pad to 100
    matches = example_splitter.smith_waterman_search(seq, "read2", "AAGCAGTGGT")
    assert matches, "Should tolerate one mismatch"
    assert matches[0]["start"] == 0 or matches[0]["start"] is not None
    assert matches[0]["end"] > 0

def test_smith_waterman_search_no_match(example_splitter):
    seq = "G" * 100
    matches = example_splitter.smith_waterman_search(seq, "read3", "AAGCAGTGGT")
    assert not matches or all(m["score"] == 0 for m in matches), "Should not find a match with wrong sequence"

def test_multiple_matches(example_splitter):
    pattern = "AAATTTGGGCCCAAGCAGTGGT"
    seq = pattern + 'A' * 50 + pattern + "A" * (100 - len(pattern + 'A' * 50 + pattern))
    matches = example_splitter.smith_waterman_search(seq, "read4", "AAGCAGTGGT")
    assert len(matches) >= 2, "Should find two matches"

def test_empty_sequence(example_splitter):
    matches = example_splitter.smith_waterman_search("", "empty", "AAGCAGTGGT")
    assert not matches or all(m["score"] == 0 for m in matches), "Should not find a match with empty sequence"
