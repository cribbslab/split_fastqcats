import unittest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import tempfile
import gzip
import os
from split_fastqcats import FastqSplitter

class TestFastqSplitter(unittest.TestCase):
    def setUp(self):
        # Create temporary files for testing
        self.temp_dir = tempfile.mkdtemp()
        self.input_fastq = os.path.join(self.temp_dir, "test_input.fastq.gz")
        self.processed_output = os.path.join(self.temp_dir, "processed.fastq.gz")
        self.lowqual_output = os.path.join(self.temp_dir, "lowqual.fastq.gz")
        self.bin_output = os.path.join(self.temp_dir, "bin.fastq.gz")
        self.stats_output = os.path.join(self.temp_dir, "stats.csv")
        
        # Test sequence with known structure
        self.test_seq = "AAGCAGTGGTATCAACGCAGAGTGAATCGTACGTACGTACGTACGTTTTTTTTTTTTCACTCTGCGTTGATACCACTGCTT"
        self.test_qual = [40] * len(self.test_seq)  # Dummy quality scores
        
        # Create test FASTQ file
        record = SeqRecord(
            Seq(self.test_seq),
            id="test_read",
            description="",
            letter_annotations={"phred_quality": self.test_qual}
        )
        
        with gzip.open(self.input_fastq, 'wt') as handle:
            SeqIO.write([record], handle, "fastq")
            
        # Initialize FastqSplitter
        self.forward_primer = "AAGCAGTGGTATCAACGCAGAGT"
        self.reverse_primer = "ACTCTGCGTTGATACCACTGCTT"
        self.index_dict = {'1': ['AAATTTGGGCCC', 'GGGCCCAAATTT']}
        self.splitter = FastqSplitter(self.forward_primer, self.reverse_primer, self.index_dict)

    def test_smith_waterman_search(self):
        """Test Smith-Waterman search function"""
        matches = self.splitter.smith_waterman_search(self.test_seq, self.forward_primer)
        self.assertTrue(len(matches) > 0)
        self.assertEqual(matches[0]['start'], 0)  # Should find primer at start

    def test_find_best_primer_pairs(self):
        """Test primer pair finding"""
        pairs = self.splitter.find_best_primer_pairs(self.test_seq)
        self.assertTrue(len(pairs) > 0)
        self.assertTrue('trimmed_seq' in pairs[0])

    def test_split_reads(self):
        """Test full read splitting functionality"""
        self.splitter.split_reads(
            self.input_fastq,
            self.processed_output,
            self.lowqual_output,
            self.bin_output,
            self.stats_output
        )
        
        # Check that output files were created
        self.assertTrue(os.path.exists(self.processed_output))
        self.assertTrue(os.path.exists(self.lowqual_output))
        self.assertTrue(os.path.exists(self.bin_output))
        self.assertTrue(os.path.exists(self.stats_output))

    def tearDown(self):
        # Clean up temporary files
        for file in [self.input_fastq, self.processed_output, self.lowqual_output, 
                    self.bin_output, self.stats_output]:
            if os.path.exists(file):
                os.remove(file)
        os.rmdir(self.temp_dir)

if __name__ == '__main__':
    unittest.main()
