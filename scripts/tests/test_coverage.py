import gzip
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import coverage as cov


def paf(read, ref, start, end, mapq=60, ref_len=100):
    return "\t".join(map(str, [read, 50, 0, 50, "+", ref, ref_len, start, end, 50, 50, mapq])) + "\n"


def run(lines, min_mapq=20, gz=True):
    d = tempfile.mkdtemp()
    p = os.path.join(d, "x.paf" + (".gz" if gz else ""))
    with (gzip.open if gz else open)(p, "wt") as fh:
        fh.writelines(lines)
    out = os.path.join(d, "cov.txt")
    cov.coverage(p, out, min_mapq)
    with open(out) as fh:
        return [l.rstrip("\n").split("\t") for l in fh][1:]


class CoverageTests(unittest.TestCase):

    def test_mates_with_same_name_both_count(self):
        (row,) = run([paf("r", "A", 0, 50), paf("r", "A", 50, 100)])
        self.assertEqual(row, ["A", "100", "1.00", "0.00", "100.00"])

    def test_identical_duplicate_line_counted_once(self):
        (row,) = run([paf("r", "A", 0, 50), paf("r", "A", 0, 50)])
        self.assertEqual(row[2:], ["0.50", "0.50", "50.00"])

    def test_low_mapq_excluded_but_reference_listed(self):
        (row,) = run([paf("r", "A", 0, 50, mapq=5)])
        self.assertEqual(row[2:], ["0.00", "0.00", "0.00"])

    def test_plain_file(self):
        (row,) = run([paf("r", "A", 0, 100)], gz=False)
        self.assertEqual(row[4], "100.00")


if __name__ == "__main__":
    unittest.main()
