import gzip
import io
import os
import sys
import tempfile
import unittest

sys.argv = sys.argv[:1]  # keep optparse from choking on unittest's own flags
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import renameFASTA_taxid as rft


def _tmp_file(content, suffix=".txt"):
    fh = tempfile.NamedTemporaryFile(mode="w", suffix=suffix, delete=False)
    fh.write(content)
    fh.close()
    return fh.name


class LoadTaxidDictTests(unittest.TestCase):

    def test_parses_name_taxid_pairs(self):
        fh = io.StringIO("seqA\t100\nseqB\t200\n")
        self.assertEqual(dict(rft.load_taxid_dict(fh)), {"seqA": "100", "seqB": "200"})


class LoadDataTests(unittest.TestCase):

    def test_reads_plain_file(self):
        path = _tmp_file("plain content\n")
        try:
            with rft.load_data(path) as fh:
                self.assertEqual(fh.read().strip(), "plain content")
        finally:
            os.remove(path)

    def test_reads_gzipped_file(self):
        fh = tempfile.NamedTemporaryFile(suffix=".gz", delete=False)
        fh.close()
        path = fh.name
        try:
            with gzip.open(path, "wt", encoding="latin-1") as gz:
                gz.write("gzipped content\n")
            with rft.load_data(path) as fh:
                self.assertEqual(fh.read().strip(), "gzipped content")
        finally:
            os.remove(path)

    def test_dash_returns_stdin(self):
        self.assertIs(rft.load_data("-"), sys.stdin)


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


class RunTests(unittest.TestCase):

    def setUp(self):
        self._saved = {
            "IN": rft.options.IN,
            "TX": rft.options.TX,
            "OUT": rft.options.OUT,
        }
        self._cleanup = []

    def tearDown(self):
        for key, value in self._saved.items():
            setattr(rft.options, key, value)
        for path in self._cleanup:
            if os.path.exists(path):
                os.remove(path)

    def test_exits_when_input_file_missing(self):
        rft.options.IN = "/no/such/input.fasta"
        rft.options.TX = _tmp_file("seqA\t100\n")
        self._cleanup.append(rft.options.TX)
        rft.options.OUT = os.path.join(tempfile.gettempdir(), "unused_out.fasta")

        with self.assertRaises(SystemExit) as ctx:
            rft.run()
        self.assertEqual(ctx.exception.code, 1)

    def test_exits_when_taxid_file_missing(self):
        rft.options.IN = _tmp_file(">seqA\nACGT\n", suffix=".fasta")
        self._cleanup.append(rft.options.IN)
        rft.options.TX = "/no/such/taxid.tsv"
        rft.options.OUT = os.path.join(tempfile.gettempdir(), "unused_out.fasta")

        with self.assertRaises(SystemExit) as ctx:
            rft.run()
        self.assertEqual(ctx.exception.code, 1)

    def test_end_to_end_writes_renamed_fasta(self):
        rft.options.IN = _tmp_file(">seqA desc\nACGT\n>seqB\nTTTT\n", suffix=".fasta")
        rft.options.TX = _tmp_file("seqA\t100\nseqB\t200\n")
        out_fh = tempfile.NamedTemporaryFile(suffix=".fasta", delete=False)
        out_fh.close()
        rft.options.OUT = out_fh.name
        self._cleanup.extend([rft.options.IN, rft.options.TX, rft.options.OUT])

        rft.run()

        with open(rft.options.OUT) as f:
            self.assertEqual(f.read(), ">100\nACGT\n>200\nTTTT\n")


if __name__ == "__main__":
    unittest.main()
