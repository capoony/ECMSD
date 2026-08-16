import io
import os
import sys
import unittest

sys.argv = sys.argv[:1]  # keep optparse from choking on unittest's own flags
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import renameFASTA_taxid as rft


class LoadTaxidDictTests(unittest.TestCase):

    def test_parses_name_taxid_pairs(self):
        fh = io.StringIO("seqA\t100\nseqB\t200\n")
        self.assertEqual(dict(rft.load_taxid_dict(fh)), {"seqA": "100", "seqB": "200"})


class RenameFastaTests(unittest.TestCase):

    def test_renames_header_to_taxid_and_keeps_sequence(self):
        fasta = ">seqA description here\nACGT\nACGT\n"
        out = io.StringIO()
        rft.rename_fasta(io.StringIO(fasta), out, {"seqA": "100"})
        self.assertEqual(out.getvalue(), ">100\nACGT\nACGT\n")

    def test_skips_sequence_for_name_not_in_taxid_dict(self):
        fasta = ">unknownSeq\nACGT\n>seqA\nTTTT\n"
        out = io.StringIO()
        rft.rename_fasta(io.StringIO(fasta), out, {"seqA": "100"})
        self.assertEqual(out.getvalue(), ">100\nTTTT\n")

    def test_dedups_by_taxid_keeping_only_first_sequence(self):
        # seqA and seqB both map to taxid 100 -> only the first is kept
        fasta = ">seqA\nAAAA\n>seqB\nCCCC\n"
        out = io.StringIO()
        rft.rename_fasta(io.StringIO(fasta), out, {"seqA": "100", "seqB": "100"})
        self.assertEqual(out.getvalue(), ">100\nAAAA\n")


if __name__ == "__main__":
    unittest.main()
