import gzip
import os
import sys
import tempfile
import unittest

sys.argv = sys.argv[:1]  # keep optparse from choking on unittest's own flags
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import LinkTaxonomy as lt


class GetAlignmentScoreTests(unittest.TestCase):

    def test_returns_score_when_ms_tag_present(self):
        fields = ["read1"] * 12 + ["tp:A:P", "ms:i:1234", "cg:Z:10M"]
        self.assertEqual(lt.get_alignment_score(fields), 1234)

    def test_returns_zero_when_no_ms_tag(self):
        fields = ["read1"] * 12 + ["tp:A:P", "cg:Z:10M"]
        self.assertEqual(lt.get_alignment_score(fields), 0)

    def test_returns_zero_when_no_optional_fields_at_all(self):
        fields = ["read1"] * 12
        self.assertEqual(lt.get_alignment_score(fields), 0)


class LoadDataTests(unittest.TestCase):

    def test_reads_plain_file(self):
        path = _tmp_file(["plain content"])
        try:
            with lt.load_data(path) as fh:
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
            with lt.load_data(path) as fh:
                self.assertEqual(fh.read().strip(), "gzipped content")
        finally:
            os.remove(path)

    def test_dash_returns_stdin(self):
        self.assertIs(lt.load_data("-"), sys.stdin)


def _tmp_file(lines):
    fh = tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False)
    fh.write("\n".join(lines) + "\n")
    fh.close()
    return fh.name


class LoadCoverageTests(unittest.TestCase):

    def test_filters_by_threshold_and_skips_header(self):
        path = _tmp_file([
            "reference\tref_length\tmean_coverage\tstd_coverage\tpct_covered",
            "taxA\t100\t5.0\t1.0\t80.0",
            "taxB\t100\t1.0\t0.5\t10.0",
            "taxC\t100\t9.0\t2.0\t50.0",
        ])
        try:
            passing = lt.load_coverage(path, 50)
            self.assertEqual(passing, {"taxA", "taxC"})
        finally:
            os.remove(path)

    def test_skips_malformed_lines(self):
        path = _tmp_file([
            "header",
            "taxA\t100\t5.0\t1.0\tnot-a-number",
            "taxB\t100",
            "taxC\t100\t9.0\t2.0\t99.0",
        ])
        try:
            passing = lt.load_coverage(path, 50)
            self.assertEqual(passing, {"taxC"})
        finally:
            os.remove(path)


class LoadBestHitsTests(unittest.TestCase):

    def test_filters_by_mapq_and_groups_by_read(self):
        paf = (
            "read1\t0\t0\t0\t+\trefA\t0\t0\t0\t0\t0\t30\n"
            "read1\t0\t0\t0\t+\trefB\t0\t0\t0\t0\t0\t30\n"
            "read2\t0\t0\t0\t+\trefC\t0\t0\t0\t0\t0\t5\n"
        )
        path = _tmp_file(paf.splitlines())
        try:
            winners = lt.load_best_hits(path, min_mapq=20)
            self.assertEqual(list(winners.keys()), ["read1"])
            self.assertEqual(len(winners["read1"]), 2)
        finally:
            os.remove(path)

    def test_skips_lines_with_fewer_than_12_fields(self):
        paf = (
            "read1\t0\t0\t0\t+\trefA\t0\t0\t0\t0\t0\t30\n"
            "truncated\tline\twith\tfew\tfields\n"
        )
        path = _tmp_file(paf.splitlines())
        try:
            winners = lt.load_best_hits(path, min_mapq=20)
            self.assertEqual(list(winners.keys()), ["read1"])
        finally:
            os.remove(path)


class TaxonTraceTests(unittest.TestCase):

    def setUp(self):
        lt.taxon_trace.cache_clear()
        lt.parents.clear()
        lt.rank_dict.clear()
        lt.names.clear()
        # root -> kingdom -> species, mimicking NCBI nodes.dmp/names.dmp shape
        lt.parents.update({"1": "1", "10": "1", "100": "10"})
        lt.rank_dict.update({"1": "domain", "10": "kingdom", "100": "species"})
        lt.names.update({"1": "root", "10": "Animalia", "100": "Homo sapiens"})

    def test_walks_path_from_node_to_root(self):
        node_path, name_path = lt.taxon_trace("100")
        self.assertEqual(node_path, "domain|kingdom|species")
        self.assertEqual(name_path, "root|Animalia|Homo sapiens")

    def test_exits_when_node_has_no_parent_and_is_not_root(self):
        # "999" is missing from the parents dict entirely (broken/truncated
        # nodes.dmp), so the walk to root can't continue
        lt.rank_dict["999"] = "species"
        with self.assertRaises(SystemExit):
            lt.taxon_trace("999")


class MainTaxonomyLoadTests(unittest.TestCase):
    """Covers main()'s parsing of nodes.dmp/names.dmp into the module-level
    parents/rank_dict/names dicts."""

    def setUp(self):
        lt.parents.clear()
        lt.rank_dict.clear()
        lt.names.clear()
        self._saved = {
            "Tax": lt.options.Tax,
            "Names": lt.options.Names,
            "PAF": lt.options.PAF,
            "OUT": lt.options.OUT,
            "Coverage": lt.options.Coverage,
        }

    def tearDown(self):
        for key, value in self._saved.items():
            setattr(lt.options, key, value)

    def test_parses_nodes_and_names_dmp(self):
        nodes_path = _tmp_file([
            "1\t|\t1\t|\troot\t|",
            "10\t|\t1\t|\tkingdom\t|",
            "100\t|\t10\t|\tspecies\t|",
        ])
        names_path = _tmp_file([
            "1\t|\troot\t|\t\t|\tscientific name\t|",
            "10\t|\tAnimalia\t|\t\t|\tscientific name\t|",
            "10\t|\tMetazoa\t|\t\t|\tsynonym\t|",
            "100\t|\tHomo sapiens\t|\t\t|\tscientific name\t|",
        ])
        try:
            lt.options.Tax = nodes_path
            lt.options.Names = names_path
            lt.options.PAF = "unused"
            lt.options.OUT = "unused"
            lt.options.Coverage = "unused"

            lt.main()

            self.assertEqual(lt.parents, {"1": "1", "10": "1", "100": "10"})
            self.assertEqual(
                lt.rank_dict, {"1": "root", "10": "kingdom", "100": "species"})
            # the "synonym" line for taxid 10 must not override its
            # scientific name
            self.assertEqual(
                lt.names, {"1": "root", "10": "Animalia", "100": "Homo sapiens"})
        finally:
            os.remove(nodes_path)
            os.remove(names_path)

    def test_exits_when_required_options_missing(self):
        lt.options.Tax = None
        lt.options.Names = "names.dmp"
        lt.options.PAF = "x"
        lt.options.OUT = "x"
        lt.options.Coverage = "x"
        with self.assertRaises(SystemExit):
            lt.main()


class FallbackLabelTests(unittest.TestCase):
    # name_path is ordered per lt.RANKS: domain..subspecies
    FULL = ["Eukaryota", "Animalia", "Chordata", "Mammalia", "Primates",
            "Hominidae", "Homo", "Homo sapiens", "NA"]

    def test_returns_requested_rank_when_populated(self):
        self.assertEqual(lt.fallback_label(self.FULL, "genus"), "Homo")

    def test_walks_up_to_next_populated_rank_when_na(self):
        # species unresolved (only genus known) -> falls back to genus, not "Unknown"
        path = ["Eukaryota", "Animalia", "Chordata", "Mammalia", "Primates",
                "Hominidae", "Homo", "NA", "NA"]
        self.assertEqual(lt.fallback_label(path, "species"), "Homo")

    def test_all_na_up_to_root_returns_unknown(self):
        path = ["NA"] * 9
        self.assertEqual(lt.fallback_label(path, "species"), "Unknown")

    def test_unknown_requested_rank_defaults_to_species_index(self):
        self.assertEqual(lt.fallback_label(self.FULL, "not-a-real-rank"), "Homo sapiens")


if __name__ == "__main__":
    unittest.main()
