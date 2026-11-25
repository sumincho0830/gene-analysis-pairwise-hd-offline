from __future__ import annotations

import json
from pathlib import Path
import unittest


RELEASE_ROOT = Path(__file__).resolve().parents[3]
DATA_PREFIX = "window.EXPERIMENT_DASHBOARD_DATA = "
CORRECTIONS = {
    "uncorrected_cpm",
    "combat_v1_cell_diagnosis",
    "combat_v2_cell_diagnosis_tissue",
}


def read_dashboard_data() -> dict:
    text = (RELEASE_ROOT / "dashboard-data.js").read_text(encoding="utf-8").strip()
    if not text.startswith(DATA_PREFIX) or not text.endswith(";"):
        raise ValueError("Unexpected dashboard-data.js wrapper")
    return json.loads(text[len(DATA_PREFIX) : -1])


class TranscriptomicsBatchCorrectionNavigationTests(unittest.TestCase):
    def test_top_level_batch_correction_tabs_are_present(self) -> None:
        html = (RELEASE_ROOT / "index.html").read_text(encoding="utf-8")
        self.assertIn('id="batchCorrectionTabs"', html)
        self.assertIn("Uncorrected CPM", html)
        self.assertIn("ComBat v1 · cell type + diagnosis", html)
        self.assertIn("ComBat v2 · cell type + diagnosis + tissue", html)
        self.assertIn("recordMatchesBatch", html)

    def test_all_three_versions_have_like_for_like_qc(self) -> None:
        records = [
            record
            for record in read_dashboard_data()["records"]
            if record.get("source_id") == "transcriptomics_batch_correction_qc"
            and record.get("dimensions", {}).get("strategy")
            == "stored_values"
        ]
        self.assertEqual(len(records), 234)
        counts = {
            correction: sum(
                record["dimensions"]["batch_correction"] == correction
                for record in records
            )
            for correction in CORRECTIONS
        }
        self.assertEqual(counts, {correction: 78 for correction in CORRECTIONS})
        self.assertEqual(
            {record["dimensions"]["batch_correction"] for record in records},
            CORRECTIONS,
        )
        for record in records:
            artifact_ids = {artifact["id"] for artifact in record["artifacts"]}
            self.assertTrue(
                {"before_heatmap", "before_pca", "before_pca_3d"}.issubset(
                    artifact_ids
                ),
                record["id"],
            )
            for item in [*record["artifacts"], *record["files"]]:
                self.assertTrue((RELEASE_ROOT / item["path"]).is_file(), item["path"])


if __name__ == "__main__":
    unittest.main()
