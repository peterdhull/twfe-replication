"""Numerical release checks; run with ``python -m unittest discover -s tests -v``.

These tests require only files inside this repository. Covariance matrices may
be singular: a matrix with more coordinates than draws cannot be full rank.
Matching the centered-draw Gram matrix verifies positive semidefiniteness
without imposing a positive-eigenvalue or invertibility requirement.
"""

import hashlib
import json
from pathlib import Path
import sys
import unittest

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from didfigures.screening import (annotation_summary, boolean_series,
                                  marker_status, paired_inference,
                                  screen_comparisons, sign_reversal)


def read_json(relative):
    return json.loads((ROOT / relative).read_text(encoding="utf-8-sig"))


def read_csv(relative):
    return pd.read_csv(ROOT / relative, float_precision="round_trip")


def study_id_for_doi(doi):
    rows = read_csv("data/provenance/studies.csv")
    selected = rows.loc[rows.source_url.eq("https://doi.org/" + doi), "study_id"]
    if len(selected) != 1:
        raise AssertionError((doi, selected.tolist()))
    return int(selected.iloc[0])


class ScreenRules(unittest.TestCase):
    """Small counterexamples catch boundary, sign, and identity-rule changes."""

    @classmethod
    def setUpClass(cls):
        cls.policy = read_json("config/screen.json")

    def example(self, **changes):
        row = dict(family="baseline", specification_id="S01-P001", modern_estimator="csa",
                   twfe_estimate=1., modern_estimate=-1., twfe_positive_se=1.,
                   modern_positive_se=1., se_gap=1., covariance_twfe_modern=.5,
                   twfe_G_eff=20., modern_G_eff=20., gap_G_eff=20.,
                   twfe_Lmax=.1, modern_Lmax=.1, gap_Lmax=.1,
                   gap_structural_identity=False, screen_prior_inference_warning=False)
        row.update(changes)
        return screen_comparisons(pd.DataFrame([row]), self.policy).iloc[0]

    def test_strict_hard_and_warning_boundaries(self):
        for role in ("twfe", "modern"):
            with self.subTest(role=role):
                self.assertFalse(self.example(**{role + "_G_eff": 5., role + "_Lmax": .5}).screen_exclude)
                self.assertTrue(self.example(**{role + "_G_eff": 5. - 1e-7}).screen_exclude)
                self.assertTrue(self.example(**{role + "_Lmax": .5 + 1e-7}).screen_exclude)
                self.assertFalse(self.example(**{role + "_G_eff": 5. - 5e-11}).screen_exclude)
        self.assertFalse(self.example(twfe_G_eff=10., modern_Lmax=.25, gap_G_eff=10., gap_Lmax=.25).screen_fragile)
        self.assertTrue(self.example(gap_G_eff=10. - 1e-7).screen_fragile)
        self.assertTrue(self.example(gap_Lmax=.25 + 1e-7).screen_fragile)

    def test_difference_concentration_warns_but_does_not_exclude(self):
        row = self.example(gap_G_eff=2., gap_Lmax=.8)
        self.assertFalse(row.screen_exclude)
        self.assertTrue(row.screen_fragile)
        self.assertEqual(row.fill_status, "fragile")

    def test_identity_does_not_waive_estimator_or_prior_warning(self):
        identity = dict(twfe_estimate=1., modern_estimate=1., se_gap=0., covariance_twfe_modern=1.,
                        gap_structural_identity=True, gap_G_eff=np.nan, gap_Lmax=np.nan)
        row = self.example(**identity)
        self.assertFalse(row.screen_fragile)
        self.assertEqual(row.fill_status, "hollow")
        self.assertTrue(self.example(**identity, twfe_G_eff=3.).screen_exclude)
        self.assertTrue(self.example(**identity, screen_prior_inference_warning=True).screen_fragile)

    def test_significance_is_paired_and_warning_suppressed(self):
        np.testing.assert_array_equal(marker_status(np.array([1.96, 1.96001, 100., np.nan]),
                                                   np.array([False, False, True, False])),
                                      ["hollow", "filled", "fragile", "hollow"])
        np.testing.assert_array_equal(sign_reversal(np.array([-1., 1., 0., -1.]),
                                                   np.array([1., -1., -1., 0.])),
                                      [True, True, False, False])
        row = self.example()
        self.assertEqual(row.abs_t, 2.)
        self.assertEqual(row.fill_status, "filled")
        self.assertTrue(row.significant_reversal)
        self.assertFalse(self.example(screen_prior_inference_warning=True).significant_reversal)

    def test_bad_inference_and_diagnostics_rejected(self):
        with self.assertRaisesRegex(ValueError, "Difference variance"):
            self.example(se_gap=2.)
        with self.assertRaisesRegex(ValueError, "nonnegative"):
            self.example(twfe_positive_se=-1.)
        with self.assertRaisesRegex(ValueError, "Nonfinite concentration"):
            self.example(gap_G_eff=np.nan)
        with self.assertRaisesRegex(ValueError, "invalid booleans"):
            boolean_series(pd.Series(["unknown"], name="flag"))


class FrozenRelease(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.comparisons = read_csv("data/results/comparisons.csv")
        cls.specifications = read_csv("data/results/specifications.csv")
        cls.auxiliary = read_csv("data/results/auxiliary_specifications.csv")
        cls.studies = read_csv("data/results/studies.csv")
        cls.coordinates = read_csv("data/covariance/coordinates.csv")
        cls.pre = read_csv("data/results/pooled_pre.csv")
        cls.manifest = read_json("data/covariance/blocks.json")
        cls.expected = read_json("reference/treatment_only_summary.json")
        cls.decisions = screen_comparisons(cls.comparisons, read_json("config/screen.json"))

    def test_complete_anonymous_indices_and_availability(self):
        self.assertEqual(len(self.comparisons), 1396)
        self.assertEqual(len(self.specifications), 351)
        self.assertEqual(len(self.auxiliary), 48)
        self.assertEqual(len(self.studies), 24)
        self.assertEqual(len(self.coordinates), 2957)
        self.assertEqual(len(self.pre), 393)
        self.assertEqual(int(self.pre.available.sum()), 138)
        for table, identifier in ((self.studies, "study_id"), (self.specifications, "specification_id"),
                                  (self.coordinates, "coordinate_id"), (self.comparisons, "comparison_id"),
                                  (self.pre, "diagnostic_id")):
            self.assertFalse(table[identifier].duplicated().any())
            self.assertTrue(table[identifier].notna().all())
        self.assertTrue(self.specifications.specification_id.str.fullmatch(r"S\d{2}-P\d{3}").all())
        self.assertTrue(self.auxiliary.specification_id.str.fullmatch(r"S\d{2}-A\d{3}").all())
        all_specs = set(self.specifications.specification_id) | set(self.auxiliary.specification_id)
        self.assertTrue(set(self.coordinates.specification_id) <= all_specs)
        self.assertTrue(set(self.comparisons.specification_id) <= set(self.specifications.specification_id))
        owner = pd.concat([self.specifications, self.auxiliary]).set_index("specification_id").study_id
        for table in (self.coordinates, self.comparisons, self.pre):
            np.testing.assert_array_equal(table.specification_id.map(owner), table.study_id)
        pre_available = self.pre.loc[self.pre.available]
        self.assertTrue(pre_available.coordinate_id.notna().all())
        self.assertTrue(self.pre.loc[~self.pre.available, "coordinate_id"].isna().all())
        self.assertTrue(self.pre.loc[~self.pre.available, "unavailable_reason_code"].notna().all())
        self.assertEqual(set(self.specifications.summary_period), {"static", "post"})
        self.assertTrue(self.pre.summary_period.eq("pre").all())

    def test_screen_recomputation_and_frozen_summaries(self):
        for name in ("twfe_screen_exclude", "modern_screen_exclude", "screen_exclude",
                     "screen_leverage_warning", "twfe_screen_warning", "modern_screen_warning",
                     "gap_screen_warning", "screen_fragile", "screen_status"):
            pd.testing.assert_series_equal(self.decisions[name], self.comparisons[name], check_dtype=False)
        kept = self.decisions.loc[~self.decisions.screen_exclude]
        self.assertEqual(kept.study_id.nunique(), self.expected["retained_studies"])
        self.assertEqual(kept.specification_id.nunique(), self.expected["retained_specifications"])
        summary = annotation_summary(kept).set_index(["family", "modern_estimator"])
        for family, estimators in self.expected["families"].items():
            for estimator, expected in estimators.items():
                with self.subTest(family=family, estimator=estimator):
                    group = kept.loc[kept.family.eq(family) & kept.modern_estimator.eq(estimator)]
                    actual = summary.loc[(family, estimator)]
                    for column, target in (("comparisons", "pairs"), ("studies", "studies"),
                                           ("fragile", "fragile"), ("sign_reversals", "reversals"),
                                           ("filled", "significant"), ("significant_reversals", "significant_reversals")):
                        self.assertEqual(actual[column], expected[target])
                    coordinates = group[["twfe_estimate", "modern_estimate"]].div(group.twfe_positive_se, axis=0)
                    self.assertEqual(int(coordinates.abs().gt(5.).any(axis=1).sum()), expected["twfe_se_capped"])
                    for numerator, share in (("sign_reversals", "sign_reversal_share"), ("filled", "filled_share"),
                                             ("significant_reversals", "significant_reversal_share")):
                        self.assertAlmostEqual(actual[share], actual[numerator] / len(group))

    def test_named_provenance_joins_without_changing_anonymous_ids(self):
        studies = read_csv("data/provenance/studies.csv")
        specs = read_csv("data/provenance/specifications.csv")
        labels = read_csv("data/provenance/coordinate_labels.csv")
        self.assertEqual(set(studies.study_id), set(self.studies.study_id))
        self.assertFalse(studies.study_id.duplicated().any())
        self.assertTrue(studies[["title", "journal", "publication_year", "source_url"]].notna().all().all())
        self.assertTrue(studies.source_url.str.startswith("https://doi.org/").all())
        self.assertEqual(studies.sort_values("study_id").publication_year.tolist(),
                         sorted(studies.publication_year.tolist()))
        pd.testing.assert_series_equal(studies.set_index("study_id").publication_year.sort_index(),
                                       self.studies.set_index("study_id").publication_year.sort_index())
        self.assertFalse(specs.specification_id.duplicated().any())
        self.assertTrue(specs.source_specification_id.notna().all())
        self.assertEqual(set(specs.specification_id),
                         set(self.specifications.specification_id) | set(self.auxiliary.specification_id))
        self.assertEqual(set(specs.loc[specs.role.eq("active_comparison"), "specification_id"]),
                         set(self.specifications.specification_id))
        owner = pd.concat([self.specifications, self.auxiliary]).set_index("specification_id").study_id
        np.testing.assert_array_equal(specs.specification_id.map(owner), specs.study_id)
        self.assertFalse(labels.coordinate_id.duplicated().any())
        self.assertEqual(set(labels.coordinate_id), set(self.coordinates.coordinate_id))
        self.assertTrue(labels.source_coordinate_label.notna().all())

    def test_new_studies_have_recorded_descriptions_and_static_targets(self):
        studies = read_csv("data/provenance/studies.csv")
        specs = read_csv("data/provenance/specifications.csv")
        for doi, year, count, cluster in (("10.1093/qje/qjab019", 2021, 14, "State"),
                                          ("10.1093/restud/rdab040", 2022, 2, "Municipality")):
            study = study_id_for_doi(doi)
            self.assertEqual(int(studies.loc[studies.study_id.eq(study), "publication_year"].iloc[0]), year)
            numerical = self.specifications.loc[self.specifications.study_id.eq(study)]
            self.assertEqual(len(numerical), count)
            self.assertTrue(numerical.summary_period.eq("static").all())
            descriptive = specs.loc[specs.study_id.eq(study) & specs.role.eq("active_comparison")]
            self.assertEqual(len(descriptive), count)
            self.assertTrue(descriptive[["published_source_location", "outcome_variable", "treatment_label",
                                         "unit_fixed_effects", "time_fixed_effects", "published_clustering"]].notna().all().all())
            self.assertTrue(descriptive.published_clustering.str.casefold().str.contains(cluster.casefold()).all())
            pre = self.pre.loc[self.pre.study_id.eq(study)]
            self.assertEqual(len(pre), 8 * count)
            self.assertFalse(pre.available.any())

    def test_every_covariance_block_and_pair_mapping(self):
        self.assertEqual(len(self.manifest["blocks"]), 25)
        self.assertEqual({x["study_id"] for x in self.manifest["blocks"]}, set(self.studies.study_id))
        self.assertIn("unavailable; not zero", self.manifest["cross_block_covariance"])
        self.assertIn("only within one block", self.manifest["cross_block_covariance"])
        block_ids = {x["covariance_block_id"] for x in self.manifest["blocks"]}
        self.assertEqual(set(self.coordinates.covariance_block_id), block_ids)
        self.assertTrue(set(self.comparisons.covariance_block_id) <= block_ids)
        self.assertEqual(sum(x["study_id"] == 2 for x in self.manifest["blocks"]), 2)
        for block in self.manifest["blocks"]:
            with self.subTest(block=block["covariance_block_id"]):
                path = ROOT / block["file"]
                self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(), block["sha256"])
                with np.load(path, allow_pickle=False) as arrays:
                    self.assertEqual(set(arrays.files), {"point", "draws", "covariance"})
                    point, draws, covariance = (arrays[x] for x in ("point", "draws", "covariance"))
                    for array in (point, draws, covariance):
                        self.assertEqual(array.dtype.kind, "f")
                        self.assertTrue(np.isfinite(array).all())
                    n = block["coordinates"]
                    self.assertEqual(point.shape, (n,))
                    self.assertEqual(draws.shape, (999, n))
                    self.assertEqual(covariance.shape, (n, n))
                    centered = draws - draws.mean(axis=0)
                    np.testing.assert_allclose(covariance, centered.T @ centered / 998, rtol=1e-9, atol=1e-12)
                    coords = self.coordinates.loc[self.coordinates.covariance_block_id.eq(block["covariance_block_id"])].sort_values("coordinate_index")
                    np.testing.assert_array_equal(coords.coordinate_index, np.arange(n))
                    np.testing.assert_allclose(coords.estimate, point, rtol=1e-13, atol=1e-14)
                    np.testing.assert_allclose(coords.positive_se, np.sqrt(np.diag(covariance)), rtol=1e-12, atol=1e-14)
                    pairs = self.comparisons.loc[self.comparisons.covariance_block_id.eq(block["covariance_block_id"])]
                    x, y = pairs.twfe_coordinate_index.to_numpy(int), pairs.modern_coordinate_index.to_numpy(int)
                    if len(pairs):
                        np.testing.assert_array_equal(coords.coordinate_id.to_numpy()[x], pairs.twfe_coordinate_id)
                        np.testing.assert_array_equal(coords.coordinate_id.to_numpy()[y], pairs.modern_coordinate_id)
                        np.testing.assert_allclose(point[x], pairs.twfe_estimate, rtol=2e-8, atol=2e-9)
                        np.testing.assert_allclose(point[y], pairs.modern_estimate, rtol=2e-8, atol=2e-9)
                        np.testing.assert_allclose(draws[:, x].std(axis=0, ddof=1), pairs.twfe_positive_se, rtol=2e-8, atol=2e-10)
                        np.testing.assert_allclose(draws[:, y].std(axis=0, ddof=1), pairs.modern_positive_se, rtol=2e-8, atol=2e-10)
                        np.testing.assert_allclose((draws[:, x] - draws[:, y]).std(axis=0, ddof=1), pairs.se_gap, rtol=2e-8, atol=2e-10)
                        np.testing.assert_allclose(covariance[x, y], pairs.covariance_twfe_modern, rtol=2e-8, atol=2e-10)
                    pre = self.pre.loc[self.pre.available & self.pre.covariance_block_id.eq(block["covariance_block_id"])]
                    pi = pre.coordinate_index.to_numpy(int)
                    if len(pre):
                        np.testing.assert_array_equal(coords.coordinate_id.to_numpy()[pi], pre.coordinate_id)
                        np.testing.assert_allclose(point[pi], pre.estimate, rtol=2e-8, atol=2e-9)
                        np.testing.assert_allclose(np.sqrt(covariance[pi, pi]), pre.positive_se, rtol=2e-8, atol=2e-9)


class PooledPreRelease(unittest.TestCase):
    """Pre coefficients are separate summaries of existing event specifications."""

    @classmethod
    def setUpClass(cls):
        cls.pre = read_csv("data/results/pre_comparisons.csv")
        cls.treatment = read_csv("data/results/comparisons.csv")
        cls.all_pairs = pd.concat([cls.treatment, cls.pre], ignore_index=True)
        cls.policy = read_json("config/screen.json")
        cls.decisions = screen_comparisons(cls.all_pairs, cls.policy)
        cls.expected = read_json("reference/expected_summary.json")
        cls.coords = read_csv("data/covariance/coordinates.csv").set_index("coordinate_id")
        cls.diagnostics = read_csv("data/results/pooled_pre.csv").set_index("diagnostic_id")
        cls.source_specifications = read_csv("data/provenance/specifications.csv").set_index("specification_id").source_specification_id
        cls.study_roles = {role: study_id_for_doi(doi) for role, doi in {
            "curriculum": "10.1086/690951", "gift": "10.1257/aer.20230008",
            "munoz": "10.1093/qje/qjad032", "fiscal": "10.3982/ECTA20612",
            "credit": "10.1086/729065"}.items()}

    def test_complete_availability_and_no_extra_specifications(self):
        availability = read_csv("data/results/pre_comparison_availability.csv")
        key = ["family", "specification_id", "modern_estimator", "summary_period"]
        self.assertEqual(len(availability), 112)
        self.assertEqual(len(self.pre), 88)
        self.assertEqual(self.pre.specification_id.nunique(), 22)
        self.assertEqual(set(self.pre.study_id), set(self.study_roles.values()))
        self.assertTrue(self.pre.summary_period.eq("pre").all())
        self.assertTrue(self.pre.specification_type.eq("event_study_pre_average").all())
        self.assertFalse(self.all_pairs.duplicated(key).any())
        self.assertEqual(self.all_pairs.specification_id.nunique(), 351)
        self.assertEqual(len(self.all_pairs), 1484)
        event_specs = set(self.treatment.loc[self.treatment.summary_period.eq("post"), "specification_id"])
        self.assertEqual(len(event_specs), 28)
        self.assertEqual(set(availability.specification_id), event_specs)
        self.assertEqual(set(self.pre.specification_id), set(availability.loc[availability.available, "specification_id"]))
        self.assertTrue(set(self.pre.specification_id) <= event_specs)
        np.testing.assert_array_equal(availability.available, availability.twfe_available & availability.modern_available)
        self.assertTrue(availability.loc[~availability.available, "unavailable_reason_code"].notna().all())
        self.assertEqual(set(availability.loc[~availability.available, "study_id"]), {
            study_id_for_doi("10.1093/qje/qjab023"), study_id_for_doi("10.3982/ECTA17951")})
        for (_, _), group in self.pre.groupby(["family", "modern_estimator"]):
            self.assertEqual(len(group), 22)
        self.assertTrue(self.pre.original_target_rows.isna().all())
        self.assertTrue(self.pre.modern_target_rows.isna().all())
        self.assertTrue(self.pre.retained_target_share.isna().all())

    def test_actual_common_windows_and_reference_interpretation(self):
        for (_, _), group in self.pre.groupby(["specification_id", "family"]):
            self.assertEqual(group.common_horizons.nunique(), 1)
        for row in self.pre.itertuples():
            registered, common = json.loads(row.registered_horizons), json.loads(row.common_horizons)
            self.assertTrue(common)
            self.assertEqual(len(common), row.retained_horizon_count)
            self.assertTrue(set(common) <= set(registered))
            self.assertEqual(common != registered, row.horizon_shortened)
            self.assertTrue(all(term == "le_m2" or (isinstance(term, int) and term < -1) for term in common))
            if row.study_id == self.study_roles["curriculum"]:
                self.assertEqual(common, ["le_m2"])
            elif row.study_id == self.study_roles["gift"]:
                self.assertEqual(common, list(range(-6, -1)))
            elif row.study_id == self.study_roles["munoz"]:
                source_id = self.source_specifications[row.specification_id]
                self.assertEqual(common, [-3, -2] if source_id == "TN-F2-ORIGINDESTYEAR" else [-5, -4, -3, -2])
            elif row.study_id == self.study_roles["fiscal"]:
                self.assertEqual(common, list(range(-10, -1)))  # Published bin codes, not calendar years.
            elif row.study_id == self.study_roles["credit"]:
                self.assertEqual(common, list(range(-20, -1)))
        np.testing.assert_array_equal(self.pre.reference_not_pure_untreated_baseline,
                                      self.pre.study_id.eq(self.study_roles["fiscal"]))
        self.assertEqual(int(self.pre.reference_not_pure_untreated_baseline.sum()), 16)

    def test_pre_pairs_use_saved_pre_coordinates_and_joint_draws(self):
        for block_id, pairs in self.pre.groupby("covariance_block_id"):
            with np.load(ROOT / "data/covariance" / f"{block_id}.npz", allow_pickle=False) as block:
                draws, point, covariance = (block[name] for name in ("draws", "point", "covariance"))
                x, y = pairs.twfe_coordinate_index.to_numpy(int), pairs.modern_coordinate_index.to_numpy(int)
                np.testing.assert_allclose(point[x], pairs.twfe_estimate, rtol=2e-8, atol=2e-9)
                np.testing.assert_allclose(point[y], pairs.modern_estimate, rtol=2e-8, atol=2e-9)
                np.testing.assert_allclose(draws[:, x].std(axis=0, ddof=1), pairs.twfe_positive_se, rtol=2e-8, atol=2e-10)
                np.testing.assert_allclose(draws[:, y].std(axis=0, ddof=1), pairs.modern_positive_se, rtol=2e-8, atol=2e-10)
                np.testing.assert_allclose((draws[:, x] - draws[:, y]).std(axis=0, ddof=1), pairs.se_gap, rtol=2e-8, atol=2e-10)
                np.testing.assert_allclose(covariance[x, y], pairs.covariance_twfe_modern, rtol=2e-8, atol=2e-10)
                for row in pairs.itertuples():
                    for role in ("twfe", "modern"):
                        coordinate = self.coords.loc[getattr(row, role + "_coordinate_id")]
                        diagnostic = self.diagnostics.loc[getattr(row, role + "_pre_diagnostic_id")]
                        self.assertEqual(coordinate.covariance_block_id, block_id)
                        self.assertEqual(coordinate.coordinate_index, getattr(row, role + "_coordinate_index"))
                        self.assertEqual(coordinate.specification_id, row.specification_id)
                        self.assertEqual(coordinate.summary_period, "pre")
                        self.assertNotEqual(coordinate.estimator, "published_native")
                        self.assertTrue(diagnostic.available)
                        self.assertEqual(diagnostic.coordinate_id, getattr(row, role + "_coordinate_id"))
                        self.assertAlmostEqual(diagnostic.estimate, getattr(row, role + "_estimate"), places=10)

    def test_pre_screen_and_ordinary_aligned_bjs_identities(self):
        actual = screen_comparisons(self.pre, self.policy)
        for name in ("screen_exclude", "screen_fragile", "screen_status", "gap_structural_identity"):
            pd.testing.assert_series_equal(actual[name], self.pre[name], check_dtype=False)
        identical = actual.loc[actual.family.eq("aligned_samples") & actual.modern_estimator.eq("bjs") & actual.study_id.ne(self.study_roles["fiscal"])]
        self.assertEqual(len(identical), 18)
        self.assertTrue(identical.gap_structural_identity.all())
        self.assertFalse(identical.gap_screen_warning.any())
        self.assertFalse(identical.fill_status.eq("filled").any())
        np.testing.assert_allclose(identical.gap, 0., atol=1e-9, rtol=0)
        np.testing.assert_allclose(identical.se_gap, 0., atol=1e-9, rtol=0)

    def test_available_event_vectors_reproduce_pre_linear_contrasts(self):
        # This block retains the complete individual coefficient vectors in
        # addition to pooled summaries, including the shortened pre window.
        labels = read_csv("data/provenance/coordinate_labels.csv").set_index("coordinate_id").source_coordinate_label
        specs = read_csv("data/provenance/specifications.csv").set_index("specification_id").source_specification_id
        lookup = dict(zip(self.coords.index.map(labels), self.coords.coordinate_index))
        munoz = self.pre.loc[self.pre.study_id.eq(self.study_roles["munoz"])]
        self.assertEqual(munoz.covariance_block_id.nunique(), 1)
        with np.load(ROOT / "data/covariance" / (munoz.covariance_block_id.iloc[0] + ".npz"), allow_pickle=False) as block:
            for row in munoz.itertuples():
                terms = json.loads(row.common_horizons)
                for role in ("twfe", "modern"):
                    estimator = row.modern_estimator if role == "modern" else "twfe"
                    descriptor = ("matched_coefficient::pre::" + row.modern_estimator
                                  if row.family == "aligned_samples" and role == "twfe"
                                  else "coefficient::" + estimator)
                    stem = specs[row.specification_id] + "::" + descriptor + "::"
                    indices = [lookup[stem + str(term)] for term in terms]
                    summary_index = getattr(row, role + "_coordinate_index")
                    np.testing.assert_allclose(block["point"][indices].mean(), block["point"][summary_index], atol=1e-12, rtol=1e-12)
                    np.testing.assert_allclose(block["draws"][:, indices].mean(axis=1), block["draws"][:, summary_index], atol=1e-12, rtol=1e-12)
                    # Pool the covariance matrix, not the marginal SEs.
                    contrast_variance = block["covariance"][np.ix_(indices, indices)].mean()
                    self.assertAlmostEqual(contrast_variance, block["covariance"][summary_index, summary_index], places=12)

    def test_combined_annotations_and_separate_summary_periods(self):
        kept = self.decisions.loc[~self.decisions.screen_exclude]
        summary = annotation_summary(kept).set_index(["family", "modern_estimator"])
        self.assertEqual(kept.specification_id.nunique(), self.expected["retained_specifications"])
        for family, estimators in self.expected["families"].items():
            for estimator, expected in estimators.items():
                group = kept.loc[kept.family.eq(family) & kept.modern_estimator.eq(estimator)]
                actual = summary.loc[(family, estimator)]
                for column, target in (("comparisons", "pairs"), ("studies", "studies"),
                                       ("fragile", "fragile"), ("sign_reversals", "reversals"),
                                       ("filled", "significant"), ("significant_reversals", "significant_reversals"),
                                       ("pre_averages", "pre_averages"), ("post_averages", "post_averages"),
                                       ("scalar_summaries", "scalar_summaries")):
                    self.assertEqual(actual[column], expected[target])
                self.assertEqual(actual.comparisons, actual.pre_averages + actual.post_averages + actual.scalar_summaries)
                coordinates = group[["twfe_estimate", "modern_estimate"]].div(group.twfe_positive_se, axis=0)
                self.assertEqual(int(coordinates.abs().gt(5.).any(axis=1).sum()), expected["twfe_se_capped"])
                for numerator, share in (("sign_reversals", "sign_reversal_share"), ("filled", "filled_share"),
                                         ("significant_reversals", "significant_reversal_share")):
                    self.assertAlmostEqual(actual[share], actual[numerator] / len(group))


if __name__ == "__main__":
    unittest.main()
