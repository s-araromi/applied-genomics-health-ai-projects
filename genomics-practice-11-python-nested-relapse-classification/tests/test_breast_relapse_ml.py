from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from breast_relapse_ml import (
    aggregate_repeated_predictions,
    bootstrap_metrics,
    construct_five_year_endpoint,
    feature_stability,
    nested_repeated_cv,
    sha256,
    validate_model_inputs,
)


class EndpointTests(unittest.TestCase):
    def test_60_month_endpoint_respects_censoring_and_late_events(self):
        clinical = pd.DataFrame({
            "follow_up_months": [12, 36, 61, 80, 60],
            "relapse": [1, 0, 0, 1, 1],
        })
        result = construct_five_year_endpoint(clinical)
        self.assertEqual(result["endpoint_eligible"].tolist(), [True, False, True, True, True])
        self.assertEqual(result.loc[[0, 2, 3, 4], "relapse_by_60_months"].tolist(), [1, 0, 0, 1])

    def test_invalid_relapse_is_rejected(self):
        with self.assertRaises(ValueError):
            construct_five_year_endpoint(pd.DataFrame({"follow_up_months": [10, 20], "relapse": [0, 2]}))


class ValidationTests(unittest.TestCase):
    def test_nonfinite_expression_is_rejected(self):
        X = pd.DataFrame({"a": [1.0, np.nan], "b": [2.0, 3.0]})
        with self.assertRaises(ValueError):
            validate_model_inputs(X, pd.Series([0, 1]), ["s1", "s2"])

    def test_duplicate_samples_are_rejected(self):
        X = pd.DataFrame({"a": [1.0, 2.0]})
        with self.assertRaises(ValueError):
            validate_model_inputs(X, pd.Series([0, 1]), ["s1", "s1"])

    def test_sha256_is_stable(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "x.txt"
            path.write_bytes(b"authentic fixture\n")
            self.assertEqual(sha256(path), "891737ada18ea29b48942742d2f4e4e156f74d8c442804e245ffadcc9599f26f")


class AuthenticDataTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        root = Path(__file__).resolve().parents[1]
        cls.clinical = pd.read_csv(root / "data" / "clinical_metadata.csv")
        cls.expression = pd.read_csv(root / "data" / "expression_genes.csv.gz", index_col=0)
        cls.probes = pd.read_csv(root / "data" / "representative_probes.csv")

    def test_authentic_subset_dimensions(self):
        self.assertEqual(self.clinical.shape[0], 286)
        self.assertEqual(self.expression.shape, (12502, 286))

    def test_public_sample_accessions_reconcile(self):
        self.assertEqual(set(self.clinical["gsm"]), set(self.expression.columns))
        self.assertTrue(self.clinical["gsm"].str.fullmatch(r"GSM\d+").all())

    def test_reported_relapse_count_is_preserved(self):
        self.assertEqual(int(self.clinical["relapse"].sum()), 107)

    def test_five_year_endpoint_counts(self):
        eligible = self.clinical[self.clinical["endpoint_eligible"]]
        self.assertEqual(len(eligible), 276)
        self.assertEqual(int(eligible["relapse_by_60_months"].sum()), 95)
        self.assertEqual(int((~self.clinical["endpoint_eligible"]).sum()), 10)

    def test_representative_probe_mapping_is_unique(self):
        self.assertEqual(len(self.probes), 12502)
        self.assertFalse(self.probes["entrez_id"].duplicated().any())
        self.assertFalse(self.probes["probe_id"].duplicated().any())

    def test_expression_values_are_log_scale_and_finite(self):
        values = self.expression.to_numpy(dtype=float)
        self.assertTrue(np.isfinite(values).all())
        self.assertGreater(values.min(), 0)
        self.assertLess(values.max(), 20)


class ModellingTests(unittest.TestCase):
    def setUp(self):
        rng = np.random.default_rng(7)
        self.y = pd.Series(np.array([0, 1] * 20))
        X = rng.normal(size=(40, 120))
        X[:, :4] += self.y.to_numpy()[:, None] * 1.2
        self.X = pd.DataFrame(X, columns=[f"GENE{i}|{1000+i}" for i in range(120)])
        self.ids = [f"S{i:03d}" for i in range(40)]

    def test_nested_cv_has_one_prediction_per_repeat(self):
        pred, folds, features = nested_repeated_cv(
            self.X, self.y, self.ids, outer_splits=4, outer_repeats=2, inner_splits=2
        )
        self.assertEqual(len(pred), 80)
        self.assertEqual(len(folds), 8)
        self.assertTrue((pred.groupby("sample_id").size() == 2).all())
        self.assertTrue(features["feature"].str.contains(r"\|").all())

    def test_nested_feature_selection_never_exceeds_training_grid(self):
        _, folds, features = nested_repeated_cv(
            self.X, self.y, self.ids, outer_splits=4, outer_repeats=1, inner_splits=2
        )
        counts = features.groupby(["repeat", "fold"]).size().to_numpy()
        self.assertTrue(np.array_equal(counts, folds["best_k"].to_numpy()))

    def test_aggregate_predictions_averages_repeats(self):
        pred = pd.DataFrame({
            "sample_id": ["A", "A", "B", "B"], "outcome": [0, 0, 1, 1],
            "predicted_probability": [0.2, 0.4, 0.7, 0.9], "repeat": [1, 2, 1, 2],
            "fold": [1, 1, 1, 1],
        })
        agg = aggregate_repeated_predictions(pred).set_index("sample_id")
        self.assertAlmostEqual(agg.loc["A", "predicted_probability"], 0.3)
        self.assertEqual(agg.loc["B", "n_predictions"], 2)

    def test_bootstrap_metrics_has_ordered_intervals(self):
        y = np.array([0] * 20 + [1] * 20)
        p = np.linspace(0.05, 0.95, 40)
        result = bootstrap_metrics(y, p, iterations=100)
        self.assertTrue((result["ci_lower"] <= result["estimate"]).all())
        self.assertTrue((result["estimate"] <= result["ci_upper"]).all())

    def test_feature_stability_uses_total_outer_models(self):
        rows = pd.DataFrame({
            "feature": ["A|1", "A|1", "B|2"], "gene_symbol": ["A", "A", "B"],
            "entrez_id": ["1", "1", "2"], "coefficient": [1.0, -0.5, -1.0],
        })
        out = feature_stability(rows, total_folds=4).set_index("feature")
        self.assertAlmostEqual(out.loc["A|1", "selection_fraction"], 0.5)
        self.assertAlmostEqual(out.loc["A|1", "sign_consistency"], 0.5)


if __name__ == "__main__":
    unittest.main()
