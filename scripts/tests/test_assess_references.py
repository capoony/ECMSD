import io
import math
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


    def test_identity_pools_discriminatory_reads_only(self):
        # A: 2 discriminatory reads, NM 0 and NM 2 over 50 bp each; 2 shared reads with NM 5
        lines = [aln("d0", "A", 0, dm=20), aln("d1", "A", 50, dm=20, NM=2)]
        for i in range(2):
            lines += [aln(f"s{i}", "A", 500, nc=2, NM=5), aln(f"s{i}", "B", 500, nc=2, NM=5)]
        refs, _ = collect(lines)
        m, _, _ = ar.judge(refs["A"], 1, 0, 0)
        self.assertAlmostEqual(m["disc_ident"], 1 - 2 / 100)
        self.assertAlmostEqual(m["disc_exact"], 0.5)
        m, _, _ = ar.judge(refs["B"], 1, 0, 0)
        self.assertEqual((m["disc_ident"], m["disc_exact"]), (None, None))

    def test_identity_skips_reads_without_nm(self):
        nm_missing = aln("x", "A", 0, dm=20).replace("NM:i:0\t", "")
        refs, _ = collect([nm_missing])
        self.assertEqual((refs["A"].disc, refs["A"].dn), (1, 0))

    def test_pcr_duplicates_counted_once(self):
        # d0 and its copy d0b: same reference, start, end and strand; d1 elsewhere
        lines = [aln("d0", "A", 0, dm=20), aln("d0b", "A", 0, dm=20), aln("d1", "A", 100, dm=20)]
        refs, stats = collect(lines)
        self.assertEqual((refs["A"].compat, refs["A"].disc, refs["A"].dups, stats["duplicates"]), (2, 2, 1, 1))
        m, _, _ = ar.judge(refs["A"], 1, 0, 0)
        self.assertAlmostEqual(m["dup_frac"], 1 / 3)
        refs, stats = collect(lines, dedup=False)
        self.assertEqual((refs["A"].disc, refs["A"].dups, stats["duplicates"]), (3, 0, 0))

    def test_shared_reads_dedup_on_all_compatible_positions(self):
        # same position on A, different position on B: not a duplicate
        lines = [aln("s0", "A", 0, nc=2), aln("s0", "B", 0, nc=2),
                 aln("s1", "A", 0, nc=2), aln("s1", "B", 300, nc=2),
                 aln("s2", "A", 0, nc=2), aln("s2", "B", 300, nc=2)]
        refs, _ = collect(lines)
        self.assertEqual((refs["A"].compat, refs["A"].dups, refs["B"].compat, refs["B"].dups), (2, 1, 2, 1))

    def test_edit_distance_decline(self):
        self.assertEqual(ar.edit_distance_decline([700, 200, 60, 20]), 1.0)
        self.assertEqual(ar.edit_distance_decline([5, 0, 0, 0]), 1.0)
        self.assertAlmostEqual(ar.edit_distance_decline([200, 400, 300, 0]), 400 / 600)  # +200, -100, -300
        self.assertIsNone(ar.edit_distance_decline([0, 0, 0, 0]))
        lines = [aln("d0", "A", 0, dm=20, NM=1), aln("d1", "A", 50, dm=20, NM=1), aln("d2", "A", 100, dm=20, NM=9)]
        refs, _ = collect(lines)
        self.assertEqual(refs["A"].dhist, [0, 2, 0, 0])  # NM above 3 is left out
        m, _, _ = ar.judge(refs["A"], 1, 0, 0)
        self.assertAlmostEqual(m["decline"], 0.5)


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


class SpatialTests(unittest.TestCase):

    def test_breadth_ratio_is_one_for_non_overlapping_reads(self):
        # 10 reads tiling 500 of 10,000 bp: depth 0.05, expected breadth 1 - exp(-0.05)
        m, _, _ = ar.judge(collect([aln(f"d{i}", "A", i * 50, dm=20, ref_len=10000) for i in range(10)])[0]["A"], 1, 0, 0)
        self.assertAlmostEqual(m["breadth_ratio"], 0.05 / (1 - math.exp(-0.05)), places=6)

    def test_stacked_reads_have_low_ratio_and_high_concentration(self):
        spread = [aln(f"s{i}", "A", i * 500, dm=20, ref_len=10000) for i in range(20)]
        stack = [aln(f"t{i}", "B", 2000 + i * 5, dm=20, ref_len=10000) for i in range(20)]
        refs, _ = collect(spread + stack)
        ms, _, _ = ar.judge(refs["A"], 1, 0, 0)
        mt, ok, why = ar.judge(refs["B"], 1, 0, 0, min_window_ratio=0.5)
        self.assertGreater(ms["breadth_ratio"], 0.95)
        self.assertLess(mt["breadth_ratio"], 0.5)
        self.assertEqual((ms["disc_windows"], ms["max_1kb"]), (1.0, 0.1))
        self.assertEqual((mt["disc_windows"], mt["max_1kb"]), (0.1, 1.0))
        self.assertFalse(ok)
        self.assertIn("clustered_discriminatory_reads", why)
        self.assertAlmostEqual(mt["window_ratio"], 0.1 / (1 - 0.9 ** 20))
        self.assertTrue(ar.judge(refs["A"], 1, 0, 0, min_window_ratio=0.5)[1])

    def test_densest_window_wraps_around_the_circular_end(self):
        reads = [aln(f"e{i}", "A", 9700, 9750, dm=20, ref_len=10000) for i in range(3)]
        reads += [aln(f"f{i}", "A", 100 + i, 150 + i, dm=20, ref_len=10000) for i in range(3)]
        reads += [aln("g", "A", 5000, dm=20, ref_len=10000)]
        m, _, _ = ar.judge(collect(reads, dedup=False)[0]["A"], 1, 0, 0)
        self.assertAlmostEqual(m["max_1kb"], 6 / 7)

    def test_no_discriminatory_reads_gives_na(self):
        m, _, _ = ar.judge(collect([aln("s", "A", 0, nc=2), aln("s", "B", 0, nc=2)])[0]["A"], 1, 0, 0)
        self.assertEqual((m["breadth_ratio"], m["disc_windows"], m["window_ratio"], m["max_1kb"]), (None,) * 4)


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
        self.assertEqual((row["7"]["discriminatory_identity"], row["7"]["discriminatory_exact"],
                          row["7"]["edit_distance_decline"]), ("1.0000", "1.000", "1.000"))
        self.assertEqual((row["8"]["discriminatory_identity"], row["8"]["edit_distance_decline"]), ("NA", "NA"))


def self_lines():
    lines = [aln(f"d{i}", "7", i * 60, dm=20) for i in range(5)]
    lines += [aln("s", "7", 700, nc=2), aln("s", "8", 700, nc=2)]
    return lines


class CallTests(unittest.TestCase):

    @staticmethod
    def calls(lines, **kw):
        refs, _ = collect(lines)
        metrics, cred = {}, []
        for r, ev in refs.items():
            m, ok, _ = ar.judge(ev, 20, 1, 0)
            metrics[r] = m
            if ok:
                cred.append(r)
        return ar.call(refs, metrics, cred, 20, **kw)

    def test_clean_spread_reads_are_confident(self):
        c = self.calls([aln(f"d{i}", "A", i * 19, dm=20) for i in range(50)])
        self.assertEqual(c["A"][1], "confident")

    def test_mismatch_hump_is_poor_fit(self):
        # most reads with 2 mismatches: not the decline expected from sequencing errors
        lines = [aln(f"d{i}", "A", i * 19, dm=20, NM=2 if i % 5 else 0) for i in range(50)]
        c = self.calls(lines)
        self.assertEqual(c["A"][1], "poor_fit")
        self.assertIn("edit_distance_decline", c["A"][3])

    def test_deep_reads_in_one_region_are_poor_fit(self):
        # 50 reads stacked on 100 bp: depth 2.5x expects ~92% breadth, observed ~10%
        c = self.calls([aln(f"d{i}", "A", i, dm=20) for i in range(50)])
        self.assertEqual(c["A"][1], "poor_fit")
        self.assertIn("breadth_ratio", c["A"][3])

    def test_references_sharing_most_reads_form_a_group(self):
        lines = []
        for i in range(100):
            lines += [aln(f"s{i}", "A", i * 9, nc=2), aln(f"s{i}", "B", i * 9, nc=2)]
        lines += [aln(f"a{i}", "A", i * 19, dm=20) for i in range(30)]
        lines += [aln(f"b{i}", "B", i * 19, dm=20) for i in range(30)]
        c = self.calls(lines)
        self.assertEqual((c["A"][1], c["B"][1]), ("ambiguous", "ambiguous"))
        self.assertEqual(c["A"][2], c["B"][2])
        self.assertTrue(c["A"][0].startswith("B:0.6"))

    def test_rare_species_sharing_conserved_reads_with_dominant_one_stays_separate(self):
        lines = [aln(f"a{i}", "A", i % 50 * 19, dm=20) for i in range(500)]
        for i in range(60):
            lines += [aln(f"s{i}", "A", i * 9, nc=2), aln(f"s{i}", "B", i * 9, nc=2)]
        lines += [aln(f"b{i}", "B", i * 19, dm=20) for i in range(25)]
        c = self.calls(lines)
        self.assertEqual(c["B"][2], "-")
        self.assertEqual(c["B"][1], "confident")

    def test_evidence_only_at_the_N_limit_is_low_evidence(self):
        c = self.calls([aln(f"d{i}", "A", i * 19, dm=20, na=21 if i < 40 else 3) for i in range(50)])
        self.assertEqual(c["A"][1], "low_evidence")


    def test_labels_are_per_taxid_unit_when_grouping(self):
        # two haplotypes of taxid 1 and a reference of taxid 2; 100 reads fit all three
        lines = []
        for i in range(100):
            lines += [aln(f"s{i}", r, i * 9, nc=3) for r in ("1|x", "1|y", "2")]
        lines += [aln(f"x{i}", "1|x", i * 19, dm=20) for i in range(30)]
        lines += [aln(f"y{i}", "1|y", i * 19 + 7, dm=20) for i in range(30)]
        lines += [aln(f"t{i}", "2", i * 19, dm=20) for i in range(30)]
        refs, _ = collect(lines, group_by_taxid=True)
        self.assertEqual(set(refs), {"1", "2"})
        metrics, cred = {}, []
        for r, ev in refs.items():
            m, ok, _ = ar.judge(ev, 20, 1, 0)
            metrics[r] = m
            cred.append(r) if ok else None
        c = ar.call(refs, metrics, cred, 20)
        # unit 1: 60 discriminatory, 160 compatible; Jaccard to unit 2 = 100 / 190 > 0.5
        self.assertEqual(refs["1"].disc, 60)
        self.assertEqual((c["1"][1], c["2"][1]), ("ambiguous", "ambiguous"))
        self.assertTrue(c["1"][0].startswith("2:0.53"))


if __name__ == "__main__":
    unittest.main()
