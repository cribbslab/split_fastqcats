import sys
from types import SimpleNamespace

# --- Mock parasail globally before any split_fastqcats import ---
class DummyResult:
    def __init__(self, score=20, end_query=10, end_ref=10):
        self.score = score
        self.end_query = end_query
        self.end_ref = end_ref
        self.ref_begin = 0
        self.query_begin = 0
        self.ref_end = end_ref
        self.query_end = end_query
        self.cigar = SimpleNamespace(string="10M")

class DummyParasail:
    @staticmethod
    def sw_scan(seq, pattern, *args, **kwargs):
        # Simulate fuzzy matching: allow up to 2 mismatches, return all matches
        max_mismatches = 2
        pattern_len = len(pattern)
        results = []
        for i in range(len(seq) - pattern_len + 1):
            window = seq[i:i+pattern_len]
            mismatches = sum(1 for a, b in zip(window, pattern) if a != b)
            if mismatches <= max_mismatches:
                results.append(DummyResult(score=pattern_len - mismatches, end_query=i+pattern_len, end_ref=pattern_len))
        return results if results else [DummyResult(score=0, end_query=0, end_ref=0)]

sys.modules["parasail"] = DummyParasail
# --- End mock ---

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
    seq = "AAATTTGGGCCCAAGCAGTGGT" + "ACGTACGTACGT" + "ACTCTGCGTT"
    matches = example_splitter.smith_waterman_search(seq, "read1")
    assert matches, "Should find a match for exact barcode+primer"
    assert matches[0]["start"] == 0 or matches[0]["start"] is not None
    assert matches[0]["end"] > 0

def test_smith_waterman_search_with_mismatch(example_splitter):
    seq = "AAATTTGGGCCAAGCAGTGGT" + "ACGTACGTACGT" + "ACTCTGCGTT"  # One C missing
    matches = example_splitter.smith_waterman_search(seq, "read2")
    assert matches, "Should tolerate one mismatch"

def test_smith_waterman_search_no_match(example_splitter):
    seq = "GGGGGGGGGGGGGGGGGGGGGGGG"
    matches = example_splitter.smith_waterman_search(seq, "read3")
    assert not matches or all(m["score"] == 0 for m in matches), "Should not find a match with wrong sequence"

def test_multiple_matches(example_splitter):
    seq = ("AAATTTGGGCCCAAGCAGTGGT" + "NNNNN" + "AAATTTGGGCCCAAGCAGTGGT")
    matches = example_splitter.smith_waterman_search(seq, "read4")
    assert len(matches) >= 2, "Should find two matches"

def test_empty_sequence(example_splitter):
    matches = example_splitter.smith_waterman_search("", "empty")
    assert matches == [] or all(m["score"] == 0 for m in matches), "Should return empty list for empty sequence"
