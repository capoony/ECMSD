"""Several references per taxid (<taxid>|<accession>) and grouping them into one unit."""
import io
import os
import sys
import tempfile
import unittest

sys.argv = sys.argv[:1]  # keep optparse from choking on unittest's own flags
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import assess_references as ar
import LinkTaxonomy as lt
import renameFASTA_taxid as rft
import select_best_alignments as sba


def paf(read, ref, start, end=None, AS=100, NM=0, tp="P", ref_len=1000, tags=()):
    end = start + 50 if end is None else end
    f = [read, 50, 0, end - start, "+", ref, ref_len, start, end, end - start, end - start, 0,
         f"NM:i:{NM}", f"AS:i:{AS}", f"tp:A:{tp}"] + list(tags)
    return "\t".join(map(str, f)) + "\n"


def tmp(text, suffix=".txt"):
    fh = tempfile.NamedTemporaryFile(mode="w", suffix=suffix, delete=False)
    fh.write(text)
    fh.close()
    return fh.name


def select(lines):
    d = tempfile.mkdtemp()
    src, dst = os.path.join(d, "in.paf"), os.path.join(d, "out.paf")
    with open(src, "w") as fh:
        fh.writelines(lines)
    sba.select_best(src, dst)
    with open(dst) as fh:
        return [l.rstrip("\n").split("\t") for l in fh]


class TaxidMarginTagTests(unittest.TestCase):

    def test_dt_is_margin_to_another_taxid_and_only_written_when_it_differs(self):
        out = select([paf("r", "7|a", 0, AS=100), paf("r", "7|b", 0, AS=98, NM=1, tp="S"),
                      paf("r", "9", 0, AS=80, NM=4, tp="S")])
        self.assertEqual(len(out), 1)
        self.assertIn("dm:i:2", out[0])
        self.assertIn("dt:i:20", out[0])
        # only haplotypes of one taxid were reported: no competitor of another taxid
        out = select([paf("s", "7|a", 0, AS=100), paf("s", "7|b", 0, AS=98, NM=1, tp="S")])
        self.assertIn("dt:i:-1", out[0])
        # one reference per taxid: dt equals dm and is not written (baseline output unchanged)
        out = select([paf("t", "7", 0, AS=100), paf("t", "9", 0, AS=80, NM=4, tp="S")])
        self.assertFalse(any(x.startswith("dt:i:") for x in out[0]))


def haplotype_reads():
    """Haplotypes 7|a (1000 bp) and 7|b (1200 bp) of one species; species 9 is lower everywhere."""
    lines = []
    for i in range(10):   # tie between the haplotypes
        lines += [paf(f"t{i}", "7|a", i * 50, tags=["nc:i:2", "dm:i:30"]),
                  paf(f"t{i}", "7|b", i * 50, ref_len=1200, tags=["nc:i:2", "dm:i:30"])]
    for i in range(10):   # best on 7|a, 2 above 7|b and 30 above species 9
        lines.append(paf(f"u{i}", "7|a", i * 50 + 3, tags=["nc:i:1", "dm:i:2", "dt:i:30"]))
    for i in range(10):   # best on 7|b, 1 above 7|a, no other taxid reported
        lines.append(paf(f"v{i}", "7|b", i * 50 + 500, ref_len=1200, tags=["nc:i:1", "dm:i:1", "dt:i:-1"]))
    return lines


class AssessGroupingTests(unittest.TestCase):

    def test_without_grouping_haplotypes_split_the_evidence(self):
        refs, _ = ar.collect(io.StringIO("".join(haplotype_reads())))
        self.assertEqual((refs["7|a"].disc, refs["7|b"].disc), (0, 0))
        self.assertEqual((refs["7|a"].weak, refs["7|b"].weak), (10, 10))

    def test_grouping_counts_one_unit_per_taxid(self):
        refs, st = ar.collect(io.StringIO("".join(haplotype_reads())), group_by_taxid=True)
        self.assertEqual(set(refs), {"7"})
        ev = refs["7"]
        self.assertEqual((ev.compat, ev.disc, ev.weak, ev.length), (30, 30, 0, 1200))
        self.assertEqual(st["shared"], 0)
        m, ok, _ = ar.judge(ev, 20, 5, 0.5)
        self.assertTrue(ok)
        self.assertAlmostEqual(m["disc_breadth"], 100.0 * 1000 / 1200)

    def test_grouping_keeps_other_taxa_apart(self):
        lines = [paf(f"s{i}", r, i * 50, tags=["nc:i:2"]) for i in range(10) for r in ("7|a", "9")]
        refs, st = ar.collect(io.StringIO("".join(lines)), group_by_taxid=True)
        self.assertEqual((refs["7"].disc, refs["9"].disc, st["shared"]), (0, 0, 10))

    def test_names_are_looked_up_by_taxid(self):
        refs, _ = ar.collect(io.StringIO("".join(haplotype_reads())))
        out = io.StringIO()
        ar.write_table(refs, out, {"7": "Drosophila simulans"}, 20, 5, 0.5)
        self.assertIn("7|a\tDrosophila_simulans", out.getvalue())


class LinkTaxonomyGroupingTests(unittest.TestCase):

    def run_link(self, group):
        d = tempfile.mkdtemp()
        nodes = tmp("1\t|\t1\t|\tno rank\t|\n5\t|\t1\t|\tgenus\t|\n7\t|\t5\t|\tspecies\t|\n9\t|\t5\t|\tspecies\t|\n")
        names = tmp("".join(f"{t}\t|\t{n}\t|\t\t|\tscientific name\t|\n"
                            for t, n in (("1", "root"), ("5", "Gen"), ("7", "Gen alpha"), ("9", "Gen beta"))))
        lines = haplotype_reads() + [paf("w0", "9", 0, tags=["nc:i:1"]), paf("w0", "7|a", 0, tags=["nc:i:2"])]
        p = tmp("".join(lines))
        cov = tmp("reference\tref_length\tmean_coverage\tstd_coverage\tpct_covered\n"
                  "7|a\t1000\t1\t0\t100\n7|b\t1200\t1\t0\t10\n9\t1000\t1\t0\t100\n")
        ev = tmp("reference\tdiscriminatory_reads\tcredible\n" +
                 ("7\t30\tyes\n9\t10\tyes\n" if group else "7|a\t0\tyes\n7|b\t0\tyes\n9\t10\tyes\n"))
        saved = dict(vars(lt.options))
        try:
            for k, v in dict(Tax=nodes, Names=names, PAF=p, Coverage=cov, Evidence=ev, CovThreshold=50,
                             MapQuality=0, TaxonomicHierarchy="species", OUT=os.path.join(d, "out"),
                             GroupByTaxid=group).items():
                setattr(lt.options, k, v)
            lt.taxon_trace.cache_clear()
            lt.run()
        finally:
            for k, v in saved.items():
                setattr(lt.options, k, v)
        with open(os.path.join(d, "out.txt")) as fh:
            rows = [l.rstrip("\n").split("\t") for l in fh][1:]
        with open(os.path.join(d, "out.ref_summary.txt")) as fh:
            summary = [l.rstrip("\n").split("\t") for l in fh][1:]
        return rows, summary

    def test_grouped_reads_go_to_the_taxid_and_shared_reads_split_by_unit(self):
        rows, summary = self.run_link(True)
        w = {}
        for r in rows:
            w[r[1]] = w.get(r[1], 0) + float(r[-1])
        # 30 haplotype reads to 7; w0 shared between 7 (30 disc) and 9 (10 disc)
        self.assertAlmostEqual(w["7"], 30.75)
        self.assertAlmostEqual(w["9"], 0.25)
        # 7|b fails the coverage threshold on its own; the unit passes through 7|a
        self.assertEqual({r[0]: r[1] for r in summary}["7|a"], "7")

    def test_ungrouped_taxid_column_is_the_taxid_of_the_reference(self):
        rows, summary = self.run_link(False)
        self.assertEqual({r[1] for r in rows}, {"7", "9"})
        # 7|b fails the coverage threshold: its 10 unique reads are unassigned
        self.assertEqual(sum(1 for r in rows if r[0].startswith("v")), 0)
        self.assertIn("7|a", {r[0] for r in summary})


class RenameAllPerTaxidTests(unittest.TestCase):

    def test_keeps_every_sequence_with_unique_names(self):
        src = ">a1 x\nAC\n>a2 y\nGT\n>b1\nTT\n>z\nCC\n"
        tx = {"a1": "7", "a2": "7", "b1": "9"}
        out = io.StringIO()
        rft.rename_fasta(io.StringIO(src), out, tx, all_per_taxid=True)
        self.assertEqual(out.getvalue(), ">7|a1\nAC\n>7|a2\nGT\n>9|b1\nTT\n")
        out = io.StringIO()
        rft.rename_fasta(io.StringIO(src), out, tx)
        self.assertEqual(out.getvalue(), ">7\nAC\n>9\nTT\n")


if __name__ == "__main__":
    unittest.main()
