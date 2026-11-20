import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import select_best_alignments as sba


def aln(read, ref, start, end, mapq=0, tp="P", AS=None, NM=None, ref_len=1000):
    f = [read, 50, 0, 50, "+", ref, ref_len, start, end, 50, end - start, mapq]
    if tp is not None:
        f.append(f"tp:A:{tp}")
    if NM is not None:
        f.append(f"NM:i:{NM}")
    if AS is not None:
        f.append(f"AS:i:{AS}")
    return "\t".join(map(str, f)) + "\n"


def run(lines, table=False, **opts):
    d = tempfile.mkdtemp()
    src, dst, tbl = (os.path.join(d, n) for n in ("in.paf", "out.paf", "t.tsv"))
    with open(src, "w") as fh:
        fh.writelines(lines)
    stats = sba.select_best(src, dst, tbl if table else None, **opts)
    with open(dst) as fh:
        out = [l.rstrip("\n").split("\t") for l in fh]
    if table:
        with open(tbl) as fh:
            return out, stats, [l.rstrip("\n").split("\t") for l in fh]
    return out, stats


class SelectBestTests(unittest.TestCase):

    def test_highest_AS_wins_even_if_secondary_and_mapq_zero(self):
        out, st = run([aln("r", "A", 0, 50, 0, "P", AS=80, NM=2),
                       aln("r", "B", 0, 50, 0, "S", AS=100, NM=0)], one_per_read=True)
        self.assertEqual([o[5] for o in out], ["B"])
        self.assertEqual(st["best_secondary"], 1)

    def test_mapq_is_not_used_for_ranking(self):
        out, _ = run([aln("r", "A", 0, 50, 60, "P", AS=70), aln("r", "B", 0, 50, 0, "S", AS=90)], one_per_read=True)
        self.assertEqual(out[0][5], "B")

    def test_AS_tie_prefers_primary_then_lower_NM(self):
        out, st = run([aln("r", "B", 0, 50, 0, "S", AS=90, NM=0), aln("r", "A", 0, 50, 0, "P", AS=90, NM=1)], one_per_read=True)
        self.assertEqual(out[0][5], "A")
        self.assertEqual(st["reads_AS_tie"], 1)
        out, _ = run([aln("r", "B", 0, 50, 0, "S", AS=90, NM=3), aln("r", "C", 0, 50, 0, "S", AS=90, NM=1)], one_per_read=True)
        self.assertEqual(out[0][5], "C")

    def test_full_tie_keeps_first_reported(self):
        out, _ = run([aln("r", "A", 0, 50, AS=90, NM=1), aln("r", "B", 0, 50, 0, "P", AS=90, NM=1)], one_per_read=True)
        self.assertEqual(out[0][5], "A")

    def test_one_row_per_read_and_mates_are_separate_reads(self):
        out, st = run([aln("p/1", "A", 0, 50, AS=100), aln("p/1", "B", 0, 50, tp="S", AS=50),
                       aln("p/2", "A", 60, 110, AS=100)], one_per_read=True)
        self.assertEqual([o[0] for o in out], ["p/1", "p/2"])
        self.assertEqual(st["reads"], 2)

    def test_missing_AS_ranks_below_present_AS(self):
        out, st = run([aln("r", "A", 0, 50, 60, "P"), aln("r", "B", 0, 50, 0, "S", AS=10)], one_per_read=True)
        self.assertEqual(out[0][5], "B")
        self.assertEqual(st["missing_AS"], 1)

    def test_read_without_any_AS_falls_back_to_primary(self):
        out, st = run([aln("r", "B", 0, 50, tp="S"), aln("r", "A", 0, 50, tp="P"),
                       aln("q", "C", 0, 50, AS=40)], one_per_read=True)
        self.assertEqual(out[0][5], "A")
        self.assertEqual(st["reads_no_AS"], 1)

    def test_missing_tp_and_NM_are_counted(self):
        _, st = run([aln("r", "A", 0, 50, tp=None, AS=40)])
        self.assertEqual((st["missing_tp"], st["missing_NM"], st["missing_AS"]), (1, 1, 0))

    def test_main_fails_when_no_AS_anywhere(self):
        d = tempfile.mkdtemp()
        src = os.path.join(d, "in.paf")
        with open(src, "w") as fh:
            fh.write(aln("r", "A", 0, 50))
        with self.assertRaises(SystemExit):
            sba.main([src, os.path.join(d, "out.paf")])

    def test_log_lines_and_blank_lines_skipped_and_all_tags_kept(self):
        out, st = run(["minimap2 log line\n", "\n", aln("r", "A", 0, 50, AS=40, NM=1)])
        self.assertEqual(st["lines"], 1)
        self.assertIn("AS:i:40", out[0])

    def test_table_lists_all_alignments_with_delta(self):
        out, _, tbl = run([aln("r", "A", 0, 50, 0, "P", AS=80, NM=2), aln("r", "B", 5, 55, 0, "S", AS=100, NM=0)],
                          table=True)
        self.assertEqual(len(tbl), 3)
        hdr = tbl[0]
        rows = {r[hdr.index("ref")]: r for r in tbl[1:]}
        self.assertEqual(rows["B"][hdr.index("best")], "1")
        self.assertEqual(rows["A"][hdr.index("AS_delta")], "20")
        self.assertEqual(rows["A"][hdr.index("n_alignments")], "2")

    # ---- compatible alignments (default mode) ----

    def test_all_tied_alignments_are_kept_with_nc_tag(self):
        out, st = run([aln("r", "A", 0, 50, AS=100), aln("r", "B", 0, 50, tp="S", AS=100),
                       aln("r", "C", 0, 50, tp="S", AS=90)])
        self.assertEqual([o[5] for o in out], ["A", "B"])
        for o in out:
            self.assertIn("nc:i:2", o)
            self.assertIn("dm:i:10", o)
        self.assertEqual((st["reads_out"], st["lines_out"], st["reads_AS_tie"]), (1, 2, 1))

    def test_unique_best_has_nc_1_and_margin_to_runner_up(self):
        out, _ = run([aln("r", "A", 0, 50, AS=100), aln("r", "B", 0, 50, tp="S", AS=94),
                      aln("r", "C", 0, 50, tp="S", AS=70)])
        self.assertEqual(len(out), 1)
        self.assertIn("nc:i:1", out[0])
        self.assertIn("dm:i:6", out[0])

    def test_no_margin_tag_without_a_lower_competitor(self):
        out, _ = run([aln("r", "A", 0, 50, AS=100)])
        self.assertIn("nc:i:1", out[0])
        self.assertFalse(any(t.startswith("dm:i:") for t in out[0]))

    def test_same_reference_twice_counts_once_for_nc(self):
        # two alignments of one read to the same reference (e.g. repeat) are one compatible reference
        out, _ = run([aln("r", "A", 0, 50, AS=100), aln("r", "A", 200, 250, tp="S", AS=100),
                      aln("r", "B", 0, 50, tp="S", AS=80)])
        self.assertEqual([o[5] for o in out], ["A"])
        self.assertIn("nc:i:1", out[0])

    def test_original_fields_unchanged_and_tags_appended(self):
        line = aln("r", "A", 0, 50, 7, "P", AS=100, NM=1)
        out, _ = run([line])
        self.assertEqual(out[0][:len(line.rstrip("\n").split("\t"))], line.rstrip("\n").split("\t"))

    # ---- quality filters ----

    def test_min_aln_len_drops_short_best_alignment(self):
        out, st = run([aln("short", "A", 0, 20, AS=40), aln("long", "A", 0, 50, AS=100)], min_aln_len=30)
        self.assertEqual([o[0] for o in out], ["long"])
        self.assertEqual(st["reads_short"], 1)

    def test_min_aln_frac_uses_query_span(self):
        # aln() gives a 50 bp read; query span is 50 here, so build a read with a short span
        f = ["r", 100, 0, 30, "+", "A", 1000, 0, 30, 30, 30, 0, "NM:i:0", "AS:i:60", "tp:A:P"]
        out, st = run(["\t".join(map(str, f)) + "\n"], min_aln_frac=0.5)
        self.assertEqual(out, [])
        self.assertEqual(st["reads_low_fraction"], 1)
        out, _ = run(["\t".join(map(str, f)) + "\n"], min_aln_frac=0.3)
        self.assertEqual(len(out), 1)

    def test_min_identity_uses_NM_over_length(self):
        out, st = run([aln("r", "A", 0, 50, AS=60, NM=10)], min_identity=0.9)
        self.assertEqual(out, [])
        self.assertEqual(st["reads_low_identity"], 1)
        out, _ = run([aln("r", "A", 0, 50, AS=60, NM=5)], min_identity=0.9)
        self.assertEqual(len(out), 1)

    def test_min_mapq_applies_to_best_mapq_among_compatible(self):
        out, st = run([aln("r", "A", 0, 50, 0, AS=100), aln("r", "B", 0, 50, 3, "S", AS=100)], min_mapq=5)
        self.assertEqual(out, [])
        self.assertEqual(st["reads_low_mapq"], 1)
        out, _ = run([aln("r", "A", 0, 50, 0, AS=100), aln("r", "B", 0, 50, 9, "S", AS=100)], min_mapq=5)
        self.assertEqual(len(out), 2)

    def test_dropped_reads_are_flagged_zero_in_table(self):
        _, _, tbl = run([aln("r", "A", 0, 20, AS=40)], table=True, min_aln_len=30)
        hdr = tbl[0]
        self.assertEqual(tbl[1][hdr.index("best")], "0")


if __name__ == "__main__":
    unittest.main()
