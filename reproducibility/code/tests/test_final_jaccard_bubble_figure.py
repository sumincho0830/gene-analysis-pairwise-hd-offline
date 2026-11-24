from __future__ import annotations

from pathlib import Path
import sys
import unittest

import numpy as np


SRC = Path(__file__).resolve().parents[1] / "src"
sys.path.insert(0, str(SRC))

from build_final_jaccard_bubble_figure import (  # noqa: E402
    benjamini_hochberg,
    observed_weighted_jaccard,
    permutation_p_values,
)


class FinalJaccardBubbleFigureTest(unittest.TestCase):
    def test_observed_matrix_uses_fuzzy_min_over_max(self) -> None:
        rna_module = np.asarray([0, 0, 1, -1], dtype=np.int32)
        rna_weight = np.asarray([0.8, 0.4, 0.6, 0.0])
        rna_sums = np.asarray([1.2, 0.6])
        protein_module = np.asarray([0, -1, 1, 0], dtype=np.int32)
        protein_weight = np.asarray([0.6, 0.0, 0.5, 0.7])
        protein_sums = np.asarray([1.3, 0.5])
        observed = observed_weighted_jaccard(
            rna_module,
            rna_weight,
            rna_sums,
            protein_module,
            protein_weight,
            protein_sums,
        )
        self.assertAlmostEqual(observed[0, 0], 0.6 / (1.2 + 1.3 - 0.6))
        self.assertAlmostEqual(observed[1, 1], 0.5 / (0.6 + 0.5 - 0.5))
        self.assertEqual(observed[0, 1], 0.0)

    def test_monte_carlo_p_values_are_deterministic_and_never_zero(self) -> None:
        rna_module = np.asarray([0, 0, 1, 1], dtype=np.int32)
        rna_weight = np.asarray([0.9, 0.8, 0.7, 0.6])
        rna_sums = np.asarray([1.7, 1.3])
        protein_module = np.asarray([0, 0, 1, 1], dtype=np.int32)
        protein_weight = np.asarray([0.9, 0.8, 0.7, 0.6])
        protein_sums = np.asarray([1.7, 1.3])
        observed = observed_weighted_jaccard(
            rna_module,
            rna_weight,
            rna_sums,
            protein_module,
            protein_weight,
            protein_sums,
        )
        first = permutation_p_values(
            observed,
            rna_module,
            rna_weight,
            rna_sums,
            protein_module,
            protein_weight,
            protein_sums,
            permutations=99,
            seed=7,
        )
        second = permutation_p_values(
            observed,
            rna_module,
            rna_weight,
            rna_sums,
            protein_module,
            protein_weight,
            protein_sums,
            permutations=99,
            seed=7,
        )
        np.testing.assert_array_equal(first, second)
        self.assertGreaterEqual(float(first.min()), 1 / 100)
        self.assertLessEqual(float(first.max()), 1.0)

    def test_bh_adjustment_matches_known_example(self) -> None:
        adjusted = benjamini_hochberg(np.asarray([0.01, 0.04, 0.03, 0.20]))
        np.testing.assert_allclose(adjusted, [0.04, 0.0533333333, 0.0533333333, 0.20])


if __name__ == "__main__":
    unittest.main()
