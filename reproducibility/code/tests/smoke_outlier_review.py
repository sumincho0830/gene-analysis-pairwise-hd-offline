"""Headless browser smoke test for the offline outlier review screen."""

from __future__ import annotations

from pathlib import Path
import tempfile

from playwright.sync_api import sync_playwright


RELEASE_ROOT = Path(__file__).resolve().parents[3]
CHROME = Path(r"C:\Program Files\Google\Chrome\Application\chrome.exe")


def main() -> None:
    console_errors: list[str] = []
    with sync_playwright() as playwright:
        browser = playwright.chromium.launch(
            executable_path=str(CHROME),
            headless=True,
        )
        page = browser.new_page(viewport={"width": 1800, "height": 1100})
        page.on(
            "console",
            lambda message: console_errors.append(message.text)
            if message.type == "error"
            else None,
        )
        page.goto((RELEASE_ROOT / "index.html").as_uri(), wait_until="load")
        page.wait_for_selector("#cellComparison .comparison-hierarchy")
        assert page.locator("#cellComparison [data-comparison-specimen]").count() == 2
        assert page.locator("#cellComparison [data-comparison-threshold]").count() == 2
        assert page.locator("#cellComparison .strategy-grid .condition-panel").count() == 3
        strategy_titles = page.locator("#cellComparison .strategy-grid .condition-head h3").all_text_contents()
        assert strategy_titles == ["CPM detection", "HVG top 9,000", "HVG top 5,000"]
        page.locator("#cellComparison [data-comparison-specimen='PB']").click()
        assert page.locator("#cellComparison [data-comparison-specimen='PB'].active").count() == 1
        page.locator(
            "#cellComparison [data-comparison-threshold='threshold_10_cells']"
        ).click()
        assert page.locator(
            "#cellComparison [data-comparison-threshold='threshold_10_cells'].active"
        ).count() == 1
        page.locator("#cellComparison [data-comparison-specimen='CSF']").click()
        page.locator(
            "#cellComparison [data-comparison-threshold='threshold_5_cells']"
        ).click()
        assert page.locator(
            "#cellComparison [data-artifact-stage='before']"
        ).count() == 1
        page.locator("#cellComparison [data-artifact-stage='before']").click()
        assert page.locator(
            "#cellComparison [data-artifact-removal]"
        ).count() == 2
        page.locator(
            "#cellComparison [data-artifact-removal='after']"
        ).click()
        page.locator(
            "#cellComparison [data-comparison-artifact='after_outlier_heatmap']"
        ).click()
        assert page.locator(
            "#cellComparison [data-comparison-artifact='after_outlier_heatmap'].active"
        ).count() == 1
        assert page.locator("#cellComparison .condition-viewer img").count() == 3

        page.locator("#cellComparison [data-artifact-stage='preprocess']").click()
        assert page.locator(
            "#cellComparison [data-artifact-removal]"
        ).count() == 2
        page.locator(
            "#cellComparison [data-artifact-removal='after']"
        ).click()
        page.locator(
            "#cellComparison "
            "[data-comparison-artifact='after_outlier_cut_diagnostic']"
        ).click()
        assert page.locator("#cellComparison .condition-viewer img").count() == 3
        page.locator(
            "#cellComparison "
            "[data-comparison-artifact='after_outlier_selected_cut']"
        ).click()
        assert page.locator("#cellComparison .condition-viewer iframe").count() == 3

        page.locator("#cellComparison [data-artifact-stage='before']").click()
        page.locator(
            "#cellComparison [data-artifact-removal='before']"
        ).click()
        page.locator(
            "#cellComparison [data-comparison-artifact='before_heatmap']"
        ).click()
        page.wait_for_selector("#cellComparison .condition-outlier-evidence")
        assert page.locator(
            "#cellComparison .condition-outlier-head", has_text="Heatmap outliers"
        ).count() == 3
        page.locator(
            "#cellComparison [data-comparison-artifact='before_pca']"
        ).click()
        page.wait_for_selector("#cellComparison .condition-outlier-evidence")
        assert page.locator("#cellComparison .condition-outlier-evidence").count() == 3
        assert page.locator(
            "#cellComparison .condition-outlier-table tbody tr"
        ).count() > 0
        evidence_text = page.locator(
            "#cellComparison .condition-outlier-evidence"
        ).first.text_content()
        assert "Distance" in evidence_text
        assert "Percentile" in evidence_text
        assert "Weighted" in evidence_text
        comparison_screenshot = (
            Path(tempfile.gettempdir()) / "gene_analysis_cell_comparison.png"
        )
        page.screenshot(path=str(comparison_screenshot), full_page=True)

        page.select_option("#viewMode", "outliers")
        page.wait_for_selector("#outlierReview .outlier-matrix-table")
        assert page.locator(
            "#outlierReview [data-outlier-dimension='stage']"
        ).count() == 5
        candidate_rows = page.locator("#outlierReview [data-outlier-sample]")
        assert candidate_rows.count() == 61
        assert page.locator("#outlierReview .outlier-cell-head").count() == 13
        assert page.locator("#outlierReview td.outlier-cell.flagged").count() > 0
        assert page.locator("#outlierReview td.outlier-cell.normal").count() > 0
        candidate_rows.first.click()
        assert page.locator("#outlierReview .outlier-cell-evidence tbody tr").count() == 13
        assert page.locator("#outlierReview .outlier-cell-evidence", has_text="Distance score").count() == 1
        assert page.locator("#outlierReview .outlier-cell-evidence", has_text="Percentile").count() == 1
        assert page.locator("#outlierReview #outlierNotes").count() == 1

        page.locator(
            "#outlierReview [data-outlier-dimension='stage']"
            "[data-outlier-value='04_variance_filtered']"
        ).click()
        page.locator(
            "#outlierReview [data-outlier-dimension='specimen']"
            "[data-outlier-value='CSF']"
        ).click()
        page.locator(
            "#outlierReview [data-outlier-dimension='threshold']"
            "[data-outlier-value='threshold_5_cells']"
        ).click()
        page.locator(
            "#outlierReview [data-outlier-dimension='strategy']"
            "[data-outlier-value='01_cpm_detection']"
        ).click()
        page.locator(
            "#outlierReview [data-outlier-dimension='phase']"
            "[data-outlier-value='02_before_preprocess']"
        ).click()
        assert page.locator(
            "#outlierReview [data-outlier-dimension].active"
        ).count() == 5
        assert page.locator("#outlierReview [data-outlier-sample]").count() == 10
        page.wait_for_selector("#outlierReview [data-open-outlier-record]")
        assert page.locator("#outlierReview [data-open-outlier-record]").count() > 0
        screenshot = Path(tempfile.gettempdir()) / "gene_analysis_outlier_review.png"
        page.screenshot(path=str(screenshot), full_page=True)
        page.locator("#outlierReview [data-open-outlier-record]").first.click()
        assert page.locator("#viewMode").input_value() == "run"
        assert page.locator("#detail .detail-head .actions", has_text="Outlier candidates").count() == 1
        page.locator("#detail [data-artifact-stage='before']").click()
        assert page.locator("#detail [data-artifact-removal]").count() == 2
        page.locator("#detail [data-artifact-removal='after']").click()
        page.locator("#detail [data-artifact='after_outlier_pca']").click()
        assert page.locator("#detail [data-artifact='after_outlier_pca'].active").count() == 1

        page.goto((RELEASE_ROOT / "index.html").as_uri(), wait_until="load")
        cut_review = page.get_by_role("button", name="Review preprocess cuts")
        assert cut_review.count() == 1
        cut_review.click()
        page.wait_for_selector("#list .item")
        assert page.locator("#list .item").count() == 168
        assert page.locator("#main .viewer img").count() == 1
        assert page.locator("#selectedCut").is_enabled()
        assert page.locator("[data-decision='use']").count() == 1
        assert page.locator("[data-decision='hold']").count() == 1
        assert page.locator("[data-decision='reject']").count() == 1
        browser.close()
    assert not console_errors, console_errors
    print(
        "OK: hierarchical cell comparison including before/after outlier "
        "removal and cut diagnostic/selected-cut states, 39 outlier experiment sets, "
        "cell-type matrix/detail, run-to-outlier navigation, and 168-record "
        "preprocess-cut review linked; "
        f"screenshots={comparison_screenshot}, {screenshot}"
    )


if __name__ == "__main__":
    main()
