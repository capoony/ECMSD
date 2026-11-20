import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import call_species as cs

CALLS = {"1": ("Drosophila mauritiana", "poor_fit", "G1", "x"),
         "2": ("Drosophila simulans", "poor_fit", "G1", "y"),
         "3": ("Drosophila melanogaster", "confident", "-", "-")}


class CallSpeciesTests(unittest.TestCase):

    def test_group_is_one_row_and_proportions_sum_to_one(self):
        rows, verdict = cs.make_calls(CALLS, {"1": 30, "2": 20, "3": 50},
                                      {"1": "Drosophila", "2": "Drosophila", "3": "Drosophila"})
        self.assertEqual(len(rows), 2)
        self.assertAlmostEqual(sum(r[2] for r in rows), 1.0)
        grp = next(r for r in rows if r[4] == "1,2")
        self.assertEqual(grp[1], "poor_fit")
        self.assertTrue(grp[0].startswith("Drosophila group (G1): D. mauritiana / D. simulans"))
        self.assertEqual(verdict, "some calls are not confident, see status")

    def test_all_poor_fit_in_one_genus_falls_back_to_genus(self):
        _, verdict = cs.make_calls(CALLS, {"1": 30, "2": 20}, {"1": "Drosophila", "2": "Drosophila"})
        self.assertIn("unsampled Drosophila species", verdict)

    def test_mixed_group_is_ambiguous(self):
        calls = dict(CALLS, **{"3": ("Drosophila melanogaster", "confident", "G1", "-")})
        rows, _ = cs.make_calls(calls, {"1": 30, "3": 50}, {"1": "Drosophila", "3": "Drosophila"})
        self.assertEqual(rows[0][1], "ambiguous")


if __name__ == "__main__":
    unittest.main()
