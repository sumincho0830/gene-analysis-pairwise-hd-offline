from __future__ import annotations

import json
from pathlib import Path
import unittest


RELEASE_ROOT = Path(__file__).resolve().parents[3]
DATA_PREFIX = "window.EXPERIMENT_DASHBOARD_DATA = "
POST_PREPROCESS_ARTIFACTS = {
    "after_outlier_after_heatmap",
    "after_outlier_after_pca",
    "after_outlier_after_pca_3d",
}


def read_dashboard_data() -> dict:
    text = (RELEASE_ROOT / "dashboard-data.js").read_text(encoding="utf-8").strip()
    if not text.startswith(DATA_PREFIX) or not text.endswith(";"):
        raise ValueError("Unexpected dashboard-data.js wrapper")
    return json.loads(text[len(DATA_PREFIX) : -1])


class AfterPreprocessOutlierNavigationTests(unittest.TestCase):
    def test_after_preprocess_stage_exposes_both_removal_states(self) -> None:
        html = (RELEASE_ROOT / "index.html").read_text(encoding="utf-8")
        self.assertIn("label:'After Preprocess'", html)
        self.assertIn("'after_heatmap','after_pca','after_pca_3d'", html)
        self.assertIn(
            "'after_outlier_after_heatmap','after_outlier_after_pca','after_outlier_after_pca_3d'",
            html,
        )
        self.assertIn("label:'Before Outlier Removal'", html)
        self.assertIn("label:'After Outlier Removal'", html)
        self.assertIn("replace(/^after_outlier_/,'')", html)

    def test_every_saved_post_preprocess_qc_set_is_packaged(self) -> None:
        data = read_dashboard_data()
        records = [
            record
            for record in data["records"]
            if record.get("source_id") == "transcriptomics_variance_findmodules"
        ]
        packaged = 0
        for record in records:
            artifact_ids = {artifact["id"] for artifact in record.get("artifacts", [])}
            has_diagnostic = "after_outlier_cut_diagnostic" in artifact_ids
            has_post_qc = POST_PREPROCESS_ARTIFACTS.issubset(artifact_ids)
            if has_diagnostic:
                self.assertTrue(has_post_qc, record["id"])
            if has_post_qc:
                packaged += 1
                file_ids = {item["id"] for item in record.get("files", [])}
                self.assertIn("after_outlier_after_pca_scores", file_ids, record["id"])
        self.assertEqual(packaged, 165)


if __name__ == "__main__":
    unittest.main()
