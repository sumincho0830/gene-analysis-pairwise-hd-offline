from __future__ import annotations

import sys
from pathlib import Path
import unittest

import numpy as np


SRC = Path(__file__).resolve().parents[1] / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from run_transcriptomics_outlier_removed_diagnostic import (  # noqa: E402
    excluded_subject,
    linkage_gap_table,
    source_sample_id,
)


class OutlierRemovedDiagnosticTests(unittest.TestCase):
    def test_confirmed_subject_matching_is_exact_and_integrated_safe(self) -> None:
        self.assertEqual(excluded_subject("PRJ-BBB-0379_Ng_MOG1_CSF"), "Ng_MOG1")
        self.assertEqual(
            excluded_subject("PRJ-BBB-0378_Ng_HD3_PB__CELL__B_cell"), "Ng_HD3"
        )
        self.assertEqual(
            source_sample_id("PRJ-BBB-0378_Ng_HD3_PB__CELL__B_cell"),
            "PRJ-BBB-0378_Ng_HD3_PB",
        )
        self.assertIsNone(excluded_subject("PRJ-BBB-0532_Heming_HD3_CSF1"))
        self.assertIsNone(excluded_subject("PRJ-BBB-0360_Ng_MOG2_CSF"))

    def test_largest_gap_suggestion_is_global_midpoint_and_review_only(self) -> None:
        sample_ids = [f"sample_{index:02d}" for index in range(5)]
        synthetic_linkage = np.asarray(
            [
                [0, 1, 1, 2],
                [2, 3, 100, 2],
                [5, 4, 101, 3],
                [6, 7, 102, 5],
            ],
            dtype=float,
        )
        table, recommendation = linkage_gap_table(
            synthetic_linkage, sample_ids
        )
        self.assertEqual(recommendation["status"], "review_required")
        self.assertEqual(recommendation["suggested_cut"], 50.5)
        self.assertEqual(recommendation["absolute_gap"], 99.0)
        self.assertEqual(
            recommendation["removed_count_if_largest_cluster_retained"], 3
        )
        self.assertTrue(recommendation["requires_caution"])
        self.assertFalse(recommendation["applied"])
        selected = table.loc[table["selected_review_suggestion"]]
        self.assertEqual(len(selected), 1)
        self.assertEqual(float(selected.iloc[0]["absolute_gap"]), 99.0)


if __name__ == "__main__":
    unittest.main()
