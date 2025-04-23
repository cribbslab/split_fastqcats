import unittest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import tempfile
import gzip
import os
from split_fastqcats import PrimerSplitter as FastqSplitter

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
        self.test_seq = "AAGCAGTGGTATCAACGCAGAGTGAATGGGCGTACGTACGTACGTACGTTTTTTTTTTTTCGTACTCTGCGTTGATACCACTGCTT"
        self.test_qual = [45] * len(self.test_seq)  # Dummy quality scores
        self.test_id = "test_read"
        
        # Create test FASTQ file
        record = SeqRecord(
            Seq(self.test_seq),
            id=self.test_id,
            description="",
            letter_annotations={"phred_quality": self.test_qual}
        )
        
        with gzip.open(self.input_fastq, 'wt') as handle:
            SeqIO.write([record], handle, "fastq")
            
        # Initialize FastqSplitter
        self.forward_primer = "AAGCAGTGGTATCAACGCAGAGTGAAT"
        self.reverse_primer = "GTACTCTGCGTTGATACCACTGCTT"
        self.error = 0.3
        self.num_workers = 4
        self.chunk_size = 1000
        self.verbose = True
        self.index_dict = {'1': ['AAATTTGGGCCC', 'GGGCCCAAATTT']} # not used for PrimerSplitter
        self.splitter = FastqSplitter(self.forward_primer, self.reverse_primer, self.error)

    def test_smith_waterman_search(self):
        """Test Smith-Waterman search function"""
        matches = self.splitter.smith_waterman_search(self.test_seq, self.test_id, self.forward_primer)
        self.assertTrue(len(matches) > 0)
        self.assertEqual(matches[0]['start'], 0)  # Should find primer at start

    def test_find_best_primer_pairs(self):
        """Test primer pair finding"""
        pairs = self.splitter.find_best_primer_pairs(self.test_seq, self.test_id)
        self.assertTrue(len(pairs) > 0)
        self.assertTrue('trimmed_seq' in pairs[0])

    def test_split_reads(self):
        """Test full read splitting functionality"""
        self.splitter.parallel_split_reads(
            self.input_fastq,
            self.processed_output,
            self.lowqual_output,
            self.bin_output,
            self.stats_output,
            self.num_workers,
            self.chunk_size
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
