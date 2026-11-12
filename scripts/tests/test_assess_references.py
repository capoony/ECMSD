import io
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import assess_references as ar

REF_LEN = 1000


def aln(read, ref, start, end=None, nc=1, dm=None, na=None, AS=100, NM=0, qlen=50, ref_len=REF_LEN):
    end = start + 50 if end is None else end
    f = [read, qlen, 0, end - start, "+", ref, ref_len, start, end, end - start, end - start, 0,
         f"NM:i:{NM}", f"AS:i:{AS}", "tp:A:P", f"nc:i:{nc}"]
    if dm is not None:
        f.append(f"dm:i:{dm}")
    if na is not None:
        f.append(f"na:i:{na}")
    return "\t".join(map(str, f)) + "\n"


def collect(lines, **kw):
    return ar.collect(io.StringIO("".join(lines)), **kw)


class CollectTests(unittest.TestCase):

    def test_shared_reads_are_compatible_but_not_discriminatory(self):
        lines = []
        for i in range(10):
            for ref in ("A", "B"):
                lines.append(aln(f"s{i}", ref, i * 50, nc=2))
        refs, st = collect(lines)
        for ref in ("A", "B"):
            self.assertEqual((refs[ref].compat, refs[ref].disc, refs[ref].weak), (10, 0, 0))
        self.assertEqual((st["shared"], st["discriminatory"]), (10, 0))

    def test_unique_best_with_margin_is_discriminatory(self):
        refs, st = collect([aln("u", "A", 0, nc=1, dm=12)])
        self.assertEqual((refs["A"].compat, refs["A"].disc), (1, 1))
        self.assertEqual(st["discriminatory"], 1)

    def test_small_margin_is_weak_not_discriminatory(self):
        refs, st = collect([aln("u", "A", 0, nc=1, dm=4)], min_margin=10)
        self.assertEqual((refs["A"].disc, refs["A"].weak), (0, 1))

    def test_no_competitor_reported_counts_as_discriminatory(self):
        refs, _ = collect([aln("u", "A", 0, nc=1)])
        self.assertEqual(refs["A"].disc, 1)

    def test_poor_alignment_quality_is_weak(self):
        # short alignment, small aligned fraction, low identity
        refs, _ = collect([aln("short", "A", 0, 20, dm=20, qlen=50)], min_aln_len=30)
        self.assertEqual((refs["A"].disc, refs["A"].weak), (0, 1))
        refs, _ = collect([aln("frac", "A", 0, 40, dm=20, qlen=100)], min_aln_frac=0.5, min_aln_len=30)
        self.assertEqual((refs["A"].disc, refs["A"].weak), (0, 1))
        refs, _ = collect([aln("ident", "A", 0, 50, dm=20, NM=10)], min_identity=0.9)
        self.assertEqual((refs["A"].disc, refs["A"].weak), (0, 1))

    def test_default_minimum_alignment_length_is_20(self):
        # a 20 bp alignment of a 20 bp read counts as evidence, a 19 bp one does not
        refs, _ = collect([aln("a", "A", 0, 20, dm=20, qlen=20), aln("b", "A", 100, 119, dm=20, qlen=19)])
        self.assertEqual((refs["A"].disc, refs["A"].weak), (1, 1))

    def test_minimum_alignment_length_can_be_changed(self):
        refs, _ = collect([aln("a", "A", 0, 20, dm=20, qlen=20)], min_aln_len=25)
        self.assertEqual((refs["A"].disc, refs["A"].weak), (0, 1))

    def test_tags_are_derived_when_absent(self):
        # PAF without nc/dm: tie and margin come from the alignments themselves
        raw = lambda r, ref, a: "\t".join(map(str, [r, 50, 0, 50, "+", ref, 1000, 0, 50, 50, 50, 0, f"AS:i:{a}"])) + "\n"
        refs, _ = collect([raw("x", "A", 100), raw("x", "B", 100), raw("y", "A", 100), raw("y", "B", 80)])
        self.assertEqual((refs["A"].compat, refs["A"].disc), (2, 1))
        self.assertEqual((refs["B"].compat, refs["B"].disc), (1, 0))

    def test_reads_at_alignment_limit_are_flagged(self):
        refs, _ = collect([aln("u", "A", 0, dm=20, na=21), aln("v", "A", 60, dm=20, na=5)], max_alns=20)
        self.assertEqual((refs["A"].disc, refs["A"].limit), (2, 1))

    def test_log_lines_are_skipped(self):
        refs, st = collect(["some minimap2 log\n", aln("u", "A", 0)])
        self.assertEqual(st["reads"], 1)


class MetricTests(unittest.TestCase):

    def test_sweep_breadth_and_depth(self):
        diff = {0: 1, 50: 1, 100: -1, 150: -1}      # [0,100) depth1 + [50,150) depth 1
        covered, total = ar.sweep(diff)
        self.assertEqual((covered, total), (150, 200))

    def test_evenness_bounds(self):
        self.assertAlmostEqual(ar.evenness([5] * 20), 1.0)
        self.assertEqual(ar.evenness([10] + [0] * 19), 0.0)
        self.assertEqual(ar.evenness([0] * 20), 0.0)
        self.assertGreaterEqual(ar.evenness([1, 0, 0]), 0.0)

    def test_breadth_of_discriminatory_reads_only(self):
        lines = [aln("d", "A", 0, 100, dm=20), aln("s", "A", 500, 600, nc=2), aln("s", "B", 500, 600, nc=2)]
        refs, _ = collect(lines)
        m, _, _ = ar.judge(refs["A"], 1, 0, 0)
        self.assertAlmostEqual(m["compat_breadth"], 20.0)
        self.assertAlmostEqual(m["disc_breadth"], 10.0)


class JudgeTests(unittest.TestCase):

    def reads_for(self, n, spread):
        step = REF_LEN // spread if spread else 0
        return [aln(f"d{i}", "A", (i % max(spread, 1)) * step, dm=20) for i in range(n)]

    def test_conserved_only_species_is_not_credible(self):
        # B: 1000 shared reads only; A: has specific reads. B has full compatible coverage but no evidence
        lines = []
        for i in range(100):
            lines += [aln(f"s{i}", "A", i * 9, nc=2), aln(f"s{i}", "B", i * 9, nc=2)]
        for i in range(40):
            lines.append(aln(f"a{i}", "A", i * 24, dm=10, AS=100))
        refs, _ = collect(lines)
        _, ok_a, _ = ar.judge(refs["A"], 20, 5, 0.5)
        _, ok_b, why_b = ar.judge(refs["B"], 20, 5, 0.5)
        self.assertTrue(ok_a)
        self.assertFalse(ok_b)
        self.assertIn("few_discriminatory_reads", why_b)

    def test_only_shared_reads_leaves_everything_unresolved(self):
        lines = []
        for i in range(100):
            lines += [aln(f"s{i}", "A", i * 9, nc=2), aln(f"s{i}", "B", i * 9, nc=2)]
        refs, _ = collect(lines)
        for ref in ("A", "B"):
            self.assertFalse(ar.judge(refs[ref], 1, 0, 0)[1])

    def test_rare_species_with_spread_evidence_is_credible(self):
        refs, _ = collect(self.reads_for(25, 25))
        m, ok, _ = ar.judge(refs["A"], 20, 5, 0.5)
        self.assertTrue(ok, m)

    def test_reads_stacked_in_one_place_are_not_credible(self):
        refs, _ = collect(self.reads_for(200, 1))
        _, ok, why = ar.judge(refs["A"], 20, 1, 0.5)
        self.assertFalse(ok)
        self.assertIn("uneven_discriminatory_coverage", why)

    def test_low_breadth_is_reported(self):
        refs, _ = collect([aln(f"d{i}", "A", 0, 10, dm=20, qlen=10) for i in range(30)])
        _, ok, why = ar.judge(refs["A"], 20, 5, 0.0)
        self.assertFalse(ok)
        self.assertIn("low_discriminatory_breadth", why)


class OutputTests(unittest.TestCase):

    def test_main_writes_table_with_names_and_flags(self):
        d = tempfile.mkdtemp()
        paf, out, names = (os.path.join(d, n) for n in ("in.paf", "out.txt", "names.dmp"))
        with open(paf, "w") as fh:
            fh.writelines(self_lines())
        with open(names, "w") as fh:
            fh.write("7\t|\tAlpha beta\t|\t\t|\tscientific name\t|\n")
        ar.main(["--PAF", paf, "--output", out, "--Names", names, "--min-disc-reads", "3",
                 "--min-disc-breadth", "1", "--min-disc-evenness", "0"])
        rows = [l.rstrip("\n").split("\t") for l in open(out)]
        hdr = rows[0]
        self.assertEqual(hdr[:2], ["reference", "taxon_name"])
        row = {r[0]: dict(zip(hdr, r)) for r in rows[1:]}
        self.assertEqual(row["7"]["taxon_name"], "Alpha_beta")
        self.assertEqual((row["7"]["credible"], row["7"]["reason"]), ("yes", "pass"))
        self.assertEqual(row["8"]["credible"], "no")
        self.assertEqual(row["8"]["taxon_name"], "NA")


def self_lines():
    lines = [aln(f"d{i}", "7", i * 60, dm=20) for i in range(5)]
    lines += [aln("s", "7", 700, nc=2), aln("s", "8", 700, nc=2)]
    return lines


if __name__ == "__main__":
    unittest.main()
