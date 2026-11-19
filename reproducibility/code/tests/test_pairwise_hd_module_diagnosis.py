from __future__ import annotations

import sys
from pathlib import Path
import unittest

import numpy as np
import pandas as pd


SRC = Path(__file__).resolve().parents[1] / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from run_pairwise_hd_module_diagnosis import calculate_pairwise_hd  # noqa: E402


class PairwiseHdModuleDiagnosisTests(unittest.TestCase):
    def test_other_diseases_are_excluded_and_direction_is_target_minus_hd(self) -> None:
        samples = [f"s{i}" for i in range(12)]
        diagnosis = pd.Series(
            ["HD"] * 3 + ["MS"] * 3 + ["MOGAD"] * 3 + ["NMOSD"] * 3,
            index=samples,
        )
        mes = pd.DataFrame(
            {
                "MEblue": [0.0, 0.1, 0.2, 1.0, 1.1, 1.2, 100, 101, 102, -100, -101, -102],
                "MEred": [1.0, 1.1, 1.2, 0.0, 0.1, 0.2, 20, 21, 22, -20, -21, -22],
            },
            index=samples,
        )
        corr, _, qval, stats, support = calculate_pairwise_hd(
            mes,
            diagnosis,
            {"blue": 30, "red": 40},
            network_id="test",
            omics="test",
            dimensions={},
        )

        self.assertGreater(corr.loc["MS vs HD", "MEblue"], 0.9)
        self.assertLess(corr.loc["MS vs HD", "MEred"], -0.9)
        self.assertTrue(np.isfinite(qval.loc["MS vs HD", "MEblue"]))
        ms = stats.loc[stats["contrast"].eq("MS vs HD")]
        self.assertTrue(ms["n_pairwise"].eq(6).all())
        self.assertEqual(int(support.loc[support["contrast"].eq("MS vs HD"), "n_pairwise"].iloc[0]), 6)

    def test_missing_group_and_constant_me_remain_na(self) -> None:
        diagnosis = pd.Series(["HD", "HD", "HD", "MS", "MS", "MS"], index=list("abcdef"))
        mes = pd.DataFrame({"MEblue": [1, 1, 1, 1, 1, 1]}, index=list("abcdef"))
        corr, pval, qval, stats, support = calculate_pairwise_hd(
            mes,
            diagnosis,
            {"blue": 6},
            network_id="test",
            omics="test",
            dimensions={},
        )

        self.assertTrue(np.isnan(corr.loc["MS vs HD", "MEblue"]))
        self.assertTrue(np.isnan(pval.loc["MOGAD vs HD", "MEblue"]))
        self.assertTrue(np.isnan(qval.loc["NMOSD vs HD", "MEblue"]))
        self.assertEqual(
            stats.loc[stats["contrast"].eq("MS vs HD"), "calculation_status"].iloc[0],
            "ZERO_VARIANCE_ME",
        )
        self.assertEqual(
            support.loc[support["contrast"].eq("MOGAD vs HD"), "support_flag"].iloc[0],
            "NO_CASE",
        )


if __name__ == "__main__":
    unittest.main()
