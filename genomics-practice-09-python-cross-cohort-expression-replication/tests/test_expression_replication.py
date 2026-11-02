"""Automated checks using disclosed genuine GEO fixture subsets."""

from __future__ import annotations

import hashlib
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT))

from expression_replication import (  # noqa: E402
    SOURCES,
    benjamini_hochberg,
    clustered_ols,
    complete_pair_indices,
    ensure_sources,
    fixed_effect_meta,
    paired_statistics,
    read_annotation,
    read_geo_matrix,
    representative_genes,
    validate_output_accounting,
)


class ExpressionReplicationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.fixtures = PROJECT / "tests/fixtures"
        cls.discovery_x, cls.discovery_m = read_geo_matrix(
            cls.fixtures / "GSE19804_four_genes.txt.gz",
            series="GSE19804", platform="GPL570", expected_shape=(4, 120),
        )
        cls.validation_x, cls.validation_m = read_geo_matrix(
            cls.fixtures / "GSE10072_four_genes.txt.gz",
            series="GSE10072", platform="GPL96", expected_shape=(4, 107),
        )

    def test_fixture_manifest_matches_authentic_subset_bytes(self) -> None:
        manifest = pd.read_csv(PROJECT / "tests/fixture_manifest.csv")
        for row in manifest.itertuples():
            path = PROJECT / "tests" / row.file
            digest = hashlib.md5(path.read_bytes(), usedforsecurity=False).hexdigest()
            self.assertEqual(path.stat().st_size, row.bytes)
            self.assertEqual(digest, row.md5)

    def test_source_manifest_is_pinned_and_unique(self) -> None:
        self.assertEqual(len(SOURCES), 4)
        self.assertEqual(len({item.name for item in SOURCES}), 4)
        for item in SOURCES:
            self.assertTrue(item.url.startswith("https://ftp.ncbi.nlm.nih.gov/geo/"))
            self.assertGreater(item.bytes, 1_000_000)
            self.assertRegex(item.md5, r"^[0-9a-f]{32}$")

    def test_missing_source_fails_when_download_disabled(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaises(FileNotFoundError):
                ensure_sources(Path(directory), allow_download=False)

    def test_authentic_cohort_dimensions_and_tissue_counts(self) -> None:
        self.assertEqual(self.discovery_x.shape, (4, 120))
        self.assertEqual(self.validation_x.shape, (4, 107))
        self.assertEqual(self.discovery_m.tissue.value_counts().to_dict(), {"Tumor": 60, "Normal": 60})
        self.assertEqual(self.validation_m.tissue.value_counts().to_dict(), {"Tumor": 58, "Normal": 49})
        self.assertEqual(self.discovery_m.patient_id.nunique(), 60)
        self.assertEqual(self.validation_m.patient_id.nunique(), 74)

    def test_complete_pairs_are_reconstructed_by_participant_id(self) -> None:
        dt, dn, dp = complete_pair_indices(self.discovery_m)
        vt, vn, vp = complete_pair_indices(self.validation_m)
        self.assertEqual((len(dp), len(vp)), (60, 33))
        self.assertTrue((self.discovery_m.loc[dt, "patient_id"].to_numpy()
                         == self.discovery_m.loc[dn, "patient_id"].to_numpy()).all())
        self.assertTrue((self.validation_m.loc[vt, "patient_id"].to_numpy()
                         == self.validation_m.loc[vn, "patient_id"].to_numpy()).all())

    def test_platform_annotation_and_gene_mapping(self) -> None:
        for platform in ("GPL570", "GPL96"):
            annotation = read_annotation(self.fixtures / f"{platform}_four_genes.annot.gz", 4)
            self.assertEqual(set(annotation.gene_symbol), {"AGER", "COL11A1", "FABP4", "SPP1"})
            expression = self.discovery_x if platform == "GPL570" else self.validation_x
            representatives = representative_genes(expression, annotation)
            self.assertEqual(len(representatives), 4)
            self.assertFalse(representatives.gene_id.duplicated().any())

    def test_paired_effects_match_genuine_fixture_counts(self) -> None:
        expected = {
            "GSE19804": np.array([-3.7473685, 4.503354, -4.397615, 3.9757295]),
            "GSE10072": np.array([-3.7975892446060606, 4.381013787757574,
                                    -4.350033962151516, 3.305101330272727]),
        }
        for series, expression, metadata in (
            ("GSE19804", self.discovery_x, self.discovery_m),
            ("GSE10072", self.validation_x, self.validation_m),
        ):
            tumor, normal, _ = complete_pair_indices(metadata)
            result = paired_statistics(expression.to_numpy(), tumor, normal)
            np.testing.assert_allclose(result.mean_log2_difference, expected[series], rtol=0, atol=1e-12)
            self.assertTrue((result.standard_error > 0).all())
            self.assertTrue((result.p_value < 1e-10).all())

    def test_pair_order_is_invariant_to_metadata_row_order(self) -> None:
        order = np.arange(len(self.validation_m))[::-1]
        shuffled_m = self.validation_m.iloc[order].reset_index(drop=True)
        shuffled_x = self.validation_x.iloc[:, order]
        tumor, normal, _ = complete_pair_indices(shuffled_m)
        result = paired_statistics(shuffled_x.to_numpy(), tumor, normal)
        original_t, original_n, _ = complete_pair_indices(self.validation_m)
        original = paired_statistics(self.validation_x.to_numpy(), original_t, original_n)
        np.testing.assert_allclose(result.mean_log2_difference, original.mean_log2_difference)

    def test_bh_values_match_independent_definition(self) -> None:
        p = np.array([0.01, 0.04, 0.03, 0.002])
        np.testing.assert_allclose(benjamini_hochberg(p), [0.02, 0.04, 0.04, 0.008])
        with self.assertRaises(ValueError):
            benjamini_hochberg(np.array([0.2, np.nan]))

    def test_clustered_ols_matches_scalar_sandwich_calculation(self) -> None:
        observed = clustered_ols(self.validation_x.to_numpy(), self.validation_m)
        x = np.column_stack((np.ones(len(self.validation_m)),
                             (self.validation_m.tissue == "Tumor").astype(float)))
        y = self.validation_x.iloc[0].to_numpy()
        inverse = np.linalg.inv(x.T @ x)
        beta = inverse @ x.T @ y
        residual = y - x @ beta
        meat = np.zeros((2, 2))
        for patient in self.validation_m.patient_id.unique():
            take = self.validation_m.patient_id.to_numpy() == patient
            score = x[take].T @ residual[take]
            meat += np.outer(score, score)
        groups, n, k = self.validation_m.patient_id.nunique(), len(y), x.shape[1]
        covariance = inverse @ meat @ inverse * groups / (groups - 1) * (n - 1) / (n - k)
        self.assertAlmostEqual(observed.clustered_effect.iloc[0], beta[1], places=12)
        self.assertAlmostEqual(observed.clustered_standard_error.iloc[0], np.sqrt(covariance[1, 1]), places=12)

    def test_fixed_effect_meta_analysis_reconciles_weights(self) -> None:
        meta = fixed_effect_meta(np.array([1.0]), np.array([0.2]), np.array([2.0]), np.array([0.5]))
        expected = (25 * 1.0 + 4 * 2.0) / 29
        self.assertAlmostEqual(meta.fixed_effect_log2_difference.iloc[0], expected)
        self.assertAlmostEqual(meta.fixed_effect_standard_error.iloc[0], np.sqrt(1 / 29))
        with self.assertRaises(ValueError):
            fixed_effect_meta(np.array([1.0]), np.array([0.0]), np.array([2.0]), np.array([0.5]))

    def test_committed_outputs_reconcile_summary_and_selection(self) -> None:
        results = pd.read_csv(PROJECT / "outputs/cross_cohort_gene_results.csv")
        replicated = pd.read_csv(PROJECT / "outputs/replicated_gene_results.csv")
        summary = pd.read_csv(PROJECT / "outputs/analysis_summary.csv").set_index("metric").value
        validate_output_accounting(results)
        self.assertEqual(len(results), int(summary.common_entrez_genes))
        self.assertEqual(int(results.discovery_candidate.sum()), int(summary.discovery_candidates))
        self.assertEqual(int(results.replicated.sum()), int(summary.replicated_genes))
        self.assertEqual(set(results.loc[results.replicated, "gene_id"]), set(replicated.gene_id))
        self.assertTrue(replicated.clustered_sensitivity_pass.all())

    def test_result_thresholds_are_applied_exactly(self) -> None:
        results = pd.read_csv(PROJECT / "outputs/cross_cohort_gene_results.csv")
        selected = results.replicated
        self.assertTrue((results.loc[selected, "discovery_q"] <= 0.05).all())
        self.assertTrue((results.loc[selected, "discovery_mean_log2_difference"].abs() >= 1.0).all())
        self.assertTrue((results.loc[selected, "validation_q_within_candidates"] <= 0.05).all())
        self.assertTrue((results.loc[selected, "validation_mean_log2_difference"].abs() >= 0.5).all())
        self.assertTrue(results.loc[selected, "validation_same_direction"].all())

    def test_wrong_series_identity_is_rejected(self) -> None:
        with self.assertRaises(ValueError):
            read_geo_matrix(self.fixtures / "GSE10072_four_genes.txt.gz",
                            series="GSE19804", platform="GPL96", expected_shape=(4, 107))


if __name__ == "__main__":
    unittest.main(verbosity=2)
