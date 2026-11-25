from __future__ import annotations

import json
from collections import Counter
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[3]


def read_assignment(path: Path) -> dict:
    text = path.read_text(encoding="utf-8")
    return json.loads(text.split("=", 1)[1].strip().rstrip(";"))


class FinalJaccardProteinScopeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.overlap = read_assignment(
            ROOT / "cross-omics-overlap/overlap-index.js"
        )["records"]
        cls.members = read_assignment(
            ROOT / "module-members/module-member-index.js"
        )["records"]
        cls.cross_index = read_assignment(
            ROOT / "module-members/cross-omics-index.js"
        )

    def test_pb_sensitivity_records_are_complete_and_linked(self) -> None:
        records = {
            record_id: entry
            for record_id, entry in self.overlap.items()
            if record_id.endswith("::protein_PB")
        }
        self.assertEqual(
            Counter(entry["batch_correction"] for entry in records.values()),
            {"combat_v1_cell_diagnosis": 138, "combat_v2_cell_diagnosis_tissue": 137},
        )
        for record_id, entry in records.items():
            self.assertEqual(entry["specimen"], "CSF")
            self.assertEqual(entry["protein_specimen"], "PB")
            self.assertNotEqual(entry["cell_type"].lower(), "other")
            self.assertIn(entry["transcript_record_id"], self.members)
            self.assertEqual(
                "::combat_v2_cell_diagnosis_tissue::" in entry["transcript_record_id"],
                entry["batch_correction"] == "combat_v2_cell_diagnosis_tissue",
            )
            self.assertTrue((ROOT / entry["path"]).is_file(), record_id)
            self.assertTrue(entry["pairs"], record_id)

    def test_pb_protein_membership_is_available_to_pair_inspector(self) -> None:
        record_id = "proteomics::PB::04_ds3_mms30_med0p10"
        modules = self.cross_index.get("proteinModules", {}).get(record_id)
        self.assertIsNotNone(modules)
        self.assertGreater(len(modules), 0)

    def test_existing_csf_records_are_explicitly_annotated(self) -> None:
        csf = [
            entry
            for record_id, entry in self.overlap.items()
            if entry.get("specimen") == "CSF"
            and not record_id.endswith("::protein_PB")
        ]
        self.assertTrue(csf)
        self.assertTrue(all(entry.get("protein_specimen") == "CSF" for entry in csf))

    def test_permutation_results_point_to_matching_version_and_diagnosis(self) -> None:
        diagnosis = read_assignment(ROOT / "module-members/cross-omics-view.js")["diagnosis"]
        for folder, version in (("", "v1"), ("combat_v2", "v2")):
            lines = (ROOT / "figures/final_jaccard_bubble" / folder / "viewer-statistics.js").read_text(
                encoding="utf-8"
            ).splitlines()
            by_threshold = json.loads(lines[-1].split("=", 1)[1].rstrip(";"))
            for threshold, statistics in by_threshold.items():
                keys = set()
                for record_id, rna, protein, p, q in statistics["pairs"]:
                    entry = self.overlap[record_id]
                    self.assertTrue(entry["batch_correction"].startswith(f"combat_{version}_"))
                    self.assertEqual(entry["threshold"], threshold)
                    self.assertEqual(entry["protein_specimen"], "CSF")
                    self.assertIn(rna, diagnosis[record_id])
                    self.assertIn(protein, entry["proteomics_modules"])
                    self.assertTrue(0 < p <= q <= 1)
                    keys.add((record_id, rna, protein))
                self.assertEqual(len(keys), statistics["module_pairs"])


if __name__ == "__main__":
    unittest.main()
