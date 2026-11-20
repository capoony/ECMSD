import gzip
import io
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import call_species as cs
import extend_refs as er
import trim_taxdump as tt

# root 1 > 7214 (genus) > 7240 (D. simulans), 7226 (D. mauritiana); 99 is unrelated
NODES = "".join(f"{c}\t|\t{p}\t|\trank\t|\n" for c, p in
                [("1", "1"), ("7214", "1"), ("7240", "7214"), ("7226", "7214"), ("99", "1")])
NAMES = "".join(f"{t}\t|\tname{t}\t|\t\t|\tscientific name\t|\n" for t in ["1", "7214", "7240", "7226", "99"])
MERGED = "5555\t|\t7226\t|\n"


def write(path, text):
    with (gzip.open(path, "wt") if path.endswith(".gz") else open(path, "w")) as fh:
        fh.write(text)
    return path


class TrimTaxdumpTests(unittest.TestCase):

    def test_keeps_lineages_of_all_fastas_and_reports_missing(self):
        with tempfile.TemporaryDirectory() as d:
            nodes = write(os.path.join(d, "nodes.dmp.gz"), NODES)
            names = write(os.path.join(d, "names.dmp.gz"), NAMES)
            ref = write(os.path.join(d, "ref.fna.gz"), ">7240\nACGT\n")
            # gzipped under a .tmp name, as extend_database() passes it
            extra = os.path.join(d, "extra.fna.gz.tmp")
            os.rename(write(os.path.join(d, "x.gz"), ">7226|AF1.1\nACGT\n>424242|X.1\nAC\n"), extra)
            out_n, out_m = os.path.join(d, "n.dmp"), os.path.join(d, "m.dmp")
            stdout = io.StringIO()
            sys.stdout, old = stdout, sys.stdout
            try:
                tt.main(["--nodes", nodes, "--names", names, "--out-nodes", out_n,
                         "--out-names", out_m, ref, extra])
            finally:
                sys.stdout = old
            kept = {l.split("\t")[0] for l in open(out_n)}
            self.assertEqual(kept, {"1", "7214", "7240", "7226"})
            self.assertEqual({l.split("\t")[0] for l in open(out_m)}, kept)
            self.assertEqual(stdout.getvalue().split(), ["424242"])
            self.assertFalse(os.path.exists(out_n + ".tmp"))


class ExtendRefsTests(unittest.TestCase):

    def setUp(self):
        self.d = tempfile.TemporaryDirectory()
        nodes = write(os.path.join(self.d.name, "nodes.dmp.gz"), NODES)
        merged = write(os.path.join(self.d.name, "merged.dmp.gz"), MERGED)
        self.known, self.remap = er.load_taxonomy(nodes, merged)

    def tearDown(self):
        self.d.cleanup()

    def run_check(self, records, skip=()):
        return [(t, i, n) for t, i, _, n in er.check(records, self.known, self.remap, set(skip))]

    def test_merged_taxid_is_replaced(self):
        self.assertEqual(self.remap, {"5555": "7226"})
        out = self.run_check([("5555|AF2.1", "A" * 6000)])
        self.assertEqual(out[0][:2], ("7226", "AF2.1"))
        self.assertIn("merged into 7226", out[0][2])

    def test_rejects_bad_header_unknown_taxid_and_duplicates(self):
        recs = [("AF3.1", "A" * 6000), ("424242|AF4.1", "A" * 6000), ("7240|AF5.1", "A" * 6000),
                ("7240|AF5.1", "A" * 6000), ("7240|OLD.1", "A" * 6000)]
        out = self.run_check(recs, skip={"OLD.1"})
        self.assertEqual([i for _, i, _ in out], ["AF5.1"])
        self.assertEqual(out[0][2], "-")

    def test_short_and_refseq_records_are_kept_with_a_note(self):
        out = self.run_check([("7240|NC_1.1", "A" * 100)])
        self.assertIn("only 100 bp", out[0][2])
        self.assertIn("RefSeq NC_", out[0][2])

    def test_read_fasta_takes_first_header_token_and_joins_lines(self):
        recs = list(er.read_fasta(["> 7240|A.1 desc", ">7240|B.1 D. simulans", "AC", "GT", ""]))
        self.assertEqual(recs, [("7240|A.1", ""), ("7240|B.1", "ACGT")])

    def test_existing_ids_strip_taxid(self):
        p = write(os.path.join(self.d.name, "extra.fna.gz"), ">7240|AF1.1\nAC\n>7226|AF2.1\nAC\n")
        self.assertEqual(er.existing_ids(p), {"AF1.1", "AF2.1"})


class CallReportTests(unittest.TestCase):

    def test_db_line_mentions_extensions_only_when_present(self):
        with tempfile.TemporaryDirectory() as d:
            p = write(os.path.join(d, "db_info.txt"),
                      "refseq_release\t238\nbuilt\t2026-11-20\ntaxdump_downloaded\t2026-11-20\n"
                      "n_extra_references\t0\nextended\tnever\n")
            self.assertEqual(cs.db_line(p), "RefSeq release 238, built 2026-11-20, taxonomy 2026-11-20")
            write(p, "refseq_release\t238\nn_extra_references\t27\nextended\t2026-11-21\n")
            self.assertIn("27 extra references (extended 2026-11-21)", cs.db_line(p))

    def test_report_warns_when_a_call_is_not_confident(self):
        rows = [["D. melanogaster", "confident", 0.6, 60, "1", "x", "-"],
                ["D. simulans", "poor_fit", 0.4, 40, "2", "y", "-"]]
        text = cs.report(rows, "some calls are not confident, see status")
        self.assertIn("poor_fit", text)
        self.assertIn("Only confident rows are species calls", text)
        self.assertNotIn("Only confident", cs.report(rows[:1], "all calls confident"))


if __name__ == "__main__":
    unittest.main()
