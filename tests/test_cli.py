"""Exercise the public entry point from an unrelated working directory."""

import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


class CommandLineReplication(unittest.TestCase):
    def test_capped_figures_and_display_coordinates(self):
        temporary_parent = ROOT / "outputs"
        temporary_parent.mkdir(exist_ok=True)
        with tempfile.TemporaryDirectory(prefix="test-", dir=temporary_parent) as temporary:
            directory = Path(temporary)
            destination = directory / "figures"
            completed = subprocess.run(
                [sys.executable, "-I", str(ROOT / "reproduce.py"), "--root", str(ROOT),
                 "--output", str(destination), "--scale", "twfe-se", "--view", "capped"],
                cwd=directory, text=True, capture_output=True, timeout=180)
            self.assertEqual(completed.returncode, 0, completed.stdout + completed.stderr)
            report = json.loads((destination / "validation.json").read_text(encoding="utf-8"))
            self.assertEqual(report["status"], "PASS")
            self.assertTrue(report["inputs_unchanged"])
            original = pd.concat([pd.read_csv(ROOT / "data/results" / name, float_precision="round_trip")
                                  for name in ("comparisons.csv", "pre_comparisons.csv")], ignore_index=True)
            coordinates = pd.read_csv(destination / "plot_coordinates.csv", float_precision="round_trip")
            joined = coordinates.merge(original, on=["family", "specification_id", "modern_estimator", "summary_period"], validate="one_to_one")
            expected = json.loads((ROOT / "reference/expected_summary.json").read_text(encoding="utf-8"))
            figure_policy = json.loads((ROOT / "config/figures.json").read_text(encoding="utf-8"))
            self.assertEqual(len(joined), sum(group["pairs"] for family in expected["families"].values() for group in family.values()))
            self.assertGreater(len(joined), 1014)
            self.assertIn("pre", set(joined.summary_period))
            self.assertEqual(len(joined.loc[~joined.summary_period.eq("pre")]), 1014)
            np.testing.assert_allclose(joined.plot_x, joined.twfe_estimate / joined.twfe_positive_se, rtol=1e-14, atol=1e-14)
            np.testing.assert_allclose(joined.plot_y, joined.modern_estimate / joined.twfe_positive_se, rtol=1e-14, atol=1e-14)
            self.assertGreater(joined.plot_x.abs().max(), 5.)  # Stored coordinates remain uncapped.
            for family in ("baseline", "aligned_samples"):
                for extension in ("png", "svg", "pdf"):
                    path = destination / "twfe-se" / f"{family}_capped.{extension}"
                    self.assertTrue(path.is_file())
                    self.assertGreater(path.stat().st_size, 1000)
                svg = (destination / "twfe-se" / f"{family}_capped.svg").read_text(encoding="utf-8")
                self.assertNotIn("CSA and BJS versus TWFE", svg)
                self.assertIn("TWFE SEs", svg)
                self.assertIn("Significant + Reversal:", svg)
                self.assertIn("Notes:", svg)
                self.assertNotIn("Variant B:", svg)
                self.assertNotIn("capped central zoom", svg)
                for panel in report["figures"]["twfe-se"][family]["capped"]:
                    self.assertEqual(panel["total"], panel["plotted"])
                    self.assertEqual(panel["total"], panel["pre_average_count"] + panel["post_average_count"] + panel["scalar_count"])
                    self.assertEqual(panel["annotation_position_axes"], figure_policy["scales"]["twfe-se"]["capped_annotation_axes"])
                    self.assertEqual(panel["annotation_point_clearance_pt"], 5)


if __name__ == "__main__":
    unittest.main()
