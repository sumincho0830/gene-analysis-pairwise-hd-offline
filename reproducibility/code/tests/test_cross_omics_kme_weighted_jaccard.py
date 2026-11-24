from __future__ import annotations

import math
from pathlib import Path
import sys
import unittest

import pandas as pd


SRC = Path(__file__).resolve().parents[1] / "src"
sys.path.insert(0, str(SRC))

from run_cross_omics_kme_weighted_jaccard import (  # noqa: E402
    compare_modules,
    module_weights,
    normalize_gene,
    weighted_jaccard,
)


class WeightedJaccardTest(unittest.TestCase):
    def test_weighted_jaccard_uses_fuzzy_union_and_intersection(self) -> None:
        weighted, plain, shared, numerator, denominator, _ = weighted_jaccard(
            {"G1": 0.8, "G2": 0.4},
            {"G1": 0.6, "G3": 0.5},
        )
        self.assertTrue(math.isclose(numerator, 0.6))
        self.assertTrue(math.isclose(denominator, 1.7))
        self.assertTrue(math.isclose(weighted, 0.6 / 1.7))
        self.assertTrue(math.isclose(plain, 1 / 3))
        self.assertEqual(shared, 1)

    def test_module_weights_keeps_modules_with_no_common_genes(self) -> None:
        table = pd.DataFrame(
            {
                "gene": ["G1", "G2"],
                "module_color": ["blue", "brown"],
                "weight": [0.8, 0.7],
                "is_raw_unassigned": [False, False],
            }
        )
        weights, counts = module_weights(
            table,
            "gene",
            "weight",
            {"G1"},
            ["blue", "brown"],
        )
        self.assertEqual(weights, {"blue": {"G1": 0.8}, "brown": {}})
        self.assertEqual(counts, {"blue": 1, "brown": 0})

    def test_gene_normalization_rejects_ambiguous_protein_symbols(self) -> None:
        self.assertEqual(normalize_gene(" hspa1a "), "HSPA1A")
        self.assertIsNone(normalize_gene("GENE1|GENE2"))
        self.assertIsNone(normalize_gene(None))

    def test_zero_common_gene_module_is_retained_but_not_called_best(self) -> None:
        pairs, _ = compare_modules(
            {"blue": {"G1": 0.8}, "snow": {}},
            {"black": {"G1": 0.7}},
            {"blue": 1, "snow": 0},
            {"black": 1},
            1,
        )
        snow = pairs.loc[pairs["transcript_module"].eq("snow")].iloc[0]
        self.assertEqual(float(snow["kme_weighted_jaccard"]), 0.0)
        self.assertFalse(bool(snow["transcript_best"]))
        self.assertFalse(bool(snow["reciprocal_best"]))

    def test_all_module_pair_shared_genes_are_returned(self) -> None:
        pairs, contributions = compare_modules(
            {"blue": {"G1": 0.8}, "red": {"G2": 0.6}},
            {"black": {"G1": 0.7}, "green": {"G2": 0.5}},
            {"blue": 1, "red": 1},
            {"black": 1, "green": 1},
            2,
        )
        self.assertEqual(len(pairs), 4)
        self.assertEqual(
            set(
                contributions[
                    ["transcript_module", "proteomics_module", "gene"]
                ].itertuples(index=False, name=None)
            ),
            {("blue", "black", "G1"), ("red", "green", "G2")},
        )


if __name__ == "__main__":
    unittest.main()
