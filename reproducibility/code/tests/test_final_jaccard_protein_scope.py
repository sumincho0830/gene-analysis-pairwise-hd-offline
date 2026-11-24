from __future__ import annotations

import json
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
        self.assertEqual(len(records), 138)
        for record_id, entry in records.items():
            self.assertEqual(entry["specimen"], "CSF")
            self.assertEqual(entry["protein_specimen"], "PB")
            self.assertNotEqual(entry["cell_type"].lower(), "other")
            self.assertIn(entry["transcript_record_id"], self.members)
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


if __name__ == "__main__":
    unittest.main()
