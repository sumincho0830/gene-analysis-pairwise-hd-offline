from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[3]


class FinalJaccardNavigationTest(unittest.TestCase):
    def test_embedded_final_heatmap_opens_parent_modal(self) -> None:
        for relative_path in (
            "final-jaccard.html",
            "reproducibility/code/src/final_jaccard_dashboard.html",
        ):
            source = (ROOT / relative_path).read_text(encoding="utf-8")
            self.assertIn("event.preventDefault()", source)
            self.assertIn("type:'open-cross-omics-overlap'", source)
            self.assertIn("document.addEventListener('click',openPairWithoutLeavingFinalPage)", source)
            self.assertIn("if(!link||window.parent===window)return", source)

    def test_parent_accepts_only_final_heatmap_messages(self) -> None:
        for relative_path in (
            "index.html",
            "reproducibility/code/src/build_experiment_inspection_dashboard.py",
        ):
            source = (ROOT / relative_path).read_text(encoding="utf-8")
            self.assertIn("event.source!==frame?.contentWindow", source)
            self.assertIn("message?.type!=='open-cross-omics-overlap'", source)
            self.assertIn("!crossEntry(message.recordId)", source)
            self.assertIn(
                "openCrossOmicsOverlap(message.recordId,message.transcriptModule||null,message.proteomicsModule||null)",
                source,
            )

    def test_final_view_can_compare_csf_and_pb_protein_networks(self) -> None:
        for relative_path in (
            "final-jaccard.html",
            "reproducibility/code/src/final_jaccard_dashboard.html",
        ):
            source = (ROOT / relative_path).read_text(encoding="utf-8")
            self.assertIn("id=\"proteinScopeControls\"", source)
            self.assertIn("['all','All protein networks']", source)
            self.assertIn("['CSF','CSF protein']", source)
            self.assertIn("['PB','PB protein']", source)
            self.assertIn("entry.protein_specimen||entry.specimen", source)

    def test_overlap_modal_uses_selected_protein_specimen(self) -> None:
        for relative_path in (
            "index.html",
            "reproducibility/code/src/build_experiment_inspection_dashboard.py",
        ):
            source = (ROOT / relative_path).read_text(encoding="utf-8")
            self.assertIn("proteinSpecimen=entry.protein_specimen||entry.specimen", source)
            self.assertIn(
                "loadModuleShard(entry.transcript_record_id||recordId,pair.transcript_module)",
                source,
            )


if __name__ == "__main__":
    unittest.main()
