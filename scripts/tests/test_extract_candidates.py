import gzip
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import extract_candidates as ec


def fq(name, seq="ACGTACGT", comment=""):
    return f"@{name}{' ' + comment if comment else ''}\n{seq}\n+\n{'I' * len(seq)}\n"


def write(path, text, gz=False):
    with (gzip.open if gz else open)(path, "wt") as fh:
        fh.write(text)


def read(path):
    with (gzip.open if path.endswith(".gz") else open)(path, "rt") as fh:
        return fh.read()


class ExtractTests(unittest.TestCase):

    def setUp(self):
        self.d = tempfile.mkdtemp()

    def p(self, name):
        return os.path.join(self.d, name)

    def run_ec(self, hits, ins, outs):
        write(self.p("hits.txt"), "".join(h + "\n" for h in hits))
        ec.main(["--hits", self.p("hits.txt"), "--in", *ins, "--out", *outs])

    def test_single_keeps_hits_unchanged_and_in_order(self):
        write(self.p("in.fq"), fq("a", "AAAA", "x y") + fq("b") + fq("c", "CCCC"))
        self.run_ec(["a", "c"], [self.p("in.fq")], [self.p("out.fq")])
        self.assertEqual(read(self.p("out.fq")), fq("a", "AAAA", "x y") + fq("c", "CCCC"))

    def test_gzip_in_and_out_detected(self):
        write(self.p("in.fq.gz"), fq("a") + fq("b"), gz=True)
        self.run_ec(["b"], [self.p("in.fq.gz")], [self.p("out.fq.gz")])
        self.assertEqual(read(self.p("out.fq.gz")), fq("b"))

    def test_gzip_input_without_gz_extension(self):
        write(self.p("in.fq"), fq("a"), gz=True)
        self.run_ec(["a"], [self.p("in.fq")], [self.p("out.fq")])
        self.assertEqual(read(self.p("out.fq")), fq("a"))

    def test_mate_suffixes_added_by_minimap2_are_ignored(self):
        # minimap2 reports 'p1/1' and 'p1/2/2' style names for mates
        write(self.p("r1.fq"), fq("p1/1") + fq("p2/1"))
        write(self.p("r2.fq"), fq("p1/2") + fq("p2/2"))
        self.run_ec(["p2/1/1"], [self.p("r1.fq"), self.p("r2.fq")], [self.p("o1.fq"), self.p("o2.fq")])
        self.assertEqual(read(self.p("o1.fq")), fq("p2/1"))
        self.assertEqual(read(self.p("o2.fq")), fq("p2/2"))

    def test_pair_kept_when_only_one_mate_hit(self):
        write(self.p("r1.fq"), fq("p1") + fq("p2"))
        write(self.p("r2.fq"), fq("p1") + fq("p2"))
        self.run_ec(["p1/2"], [self.p("r1.fq"), self.p("r2.fq")], [self.p("o1.fq"), self.p("o2.fq")])
        self.assertEqual(read(self.p("o1.fq")), fq("p1"))
        self.assertEqual(read(self.p("o2.fq")), fq("p1"))

    def test_no_hits_gives_empty_output(self):
        write(self.p("in.fq"), fq("a"))
        self.run_ec([], [self.p("in.fq")], [self.p("out.fq")])
        self.assertEqual(read(self.p("out.fq")), "")

    def test_unsynchronised_mates_rejected(self):
        write(self.p("r1.fq"), fq("p1") + fq("p2"))
        write(self.p("r2.fq"), fq("p2") + fq("p1"))
        with self.assertRaises(SystemExit):
            self.run_ec(["p1"], [self.p("r1.fq"), self.p("r2.fq")], [self.p("o1.fq"), self.p("o2.fq")])

    def test_different_read_counts_rejected(self):
        write(self.p("r1.fq"), fq("p1") + fq("p2"))
        write(self.p("r2.fq"), fq("p1"))
        with self.assertRaises(SystemExit):
            self.run_ec(["p1"], [self.p("r1.fq"), self.p("r2.fq")], [self.p("o1.fq"), self.p("o2.fq")])

    def test_truncated_fastq_rejected(self):
        write(self.p("in.fq"), fq("a") + "@b\nACGT\n")
        with self.assertRaises(SystemExit):
            self.run_ec(["a"], [self.p("in.fq")], [self.p("out.fq")])


if __name__ == "__main__":
    unittest.main()
