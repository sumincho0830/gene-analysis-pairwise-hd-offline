from __future__ import annotations

import json
from pathlib import Path
import unittest


RELEASE_ROOT = Path(__file__).resolve().parents[3]
DATA_PREFIX = "window.EXPERIMENT_DASHBOARD_DATA = "
EXPECTED_VISUALS = {
    "qc_heatmap",
    "qc_pca_metadata",
    "qc_pca_3d",
    "qc_pca_3d_interactive",
}


def read_dashboard_data() -> dict:
    text = (RELEASE_ROOT / "dashboard-data.js").read_text(encoding="utf-8").strip()
    if not text.startswith(DATA_PREFIX) or not text.endswith(";"):
        raise ValueError("Unexpected dashboard-data.js wrapper")
    return json.loads(text[len(DATA_PREFIX) : -1])


class ProteomicsQcComparisonTests(unittest.TestCase):
    def test_raw_and_batch_corrected_are_available_for_each_specimen(self) -> None:
        data = read_dashboard_data()
        records = [
            record
            for record in data["records"]
            if record.get("source_id") == "proteomics_batch_correction_qc"
        ]
        self.assertEqual(len(records), 6)
        self.assertEqual(
            {
                (record["dimensions"]["strategy"], record["dimensions"]["specimen"])
                for record in records
            },
            {
                (state, specimen)
                for state in ("00_raw", "01_batch_corrected")
                for specimen in ("ALL", "CSF", "PB")
            },
        )

        for record in records:
            artifact_ids = {artifact["id"] for artifact in record["artifacts"]}
            self.assertTrue(EXPECTED_VISUALS.issubset(artifact_ids), record["id"])
            self.assertGreater(record["metrics"]["samples"], 0, record["id"])
            self.assertGreater(record["metrics"]["genes"], 0, record["id"])
            for item in [*record["artifacts"], *record["files"]]:
                self.assertTrue((RELEASE_ROOT / item["path"]).is_file(), item["path"])

    def test_comparison_ui_orders_qc_outputs_and_labels_conditions(self) -> None:
        html = (RELEASE_ROOT / "index.html").read_text(encoding="utf-8")
        self.assertIn("'qc_heatmap','qc_pca_metadata','qc_pca_3d'", html)
        self.assertIn("'00_raw':'Raw'", html)
        self.assertIn("'01_batch_corrected':'Batch corrected'", html)


if __name__ == "__main__":
    unittest.main()
