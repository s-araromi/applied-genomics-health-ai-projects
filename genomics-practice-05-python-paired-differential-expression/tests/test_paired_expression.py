import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT))

from analyze_paired_expression import (  # noqa: E402
    ANNOTATION_BYTES,
    ANNOTATION_MD5,
    MATRIX_BYTES,
    MATRIX_MD5,
    benjamini_hochberg,
    md5sum,
    pair_indices,
    paired_test,
    unambiguous_symbol,
    validate_inputs,
)


class PairedExpressionTests(unittest.TestCase):
    def setUp(self):
        patients = ["1", "2", "3", "4"]
        self.metadata = pd.DataFrame({
            "sample_accession": [f"T{i}" for i in patients] + [f"N{i}" for i in patients],
            "sample_title": [f"Tumor {i}" for i in patients] + [f"Normal {i}" for i in patients],
            "patient_id": patients + patients,
            "tissue": ["Tumor"] * 4 + ["Normal"] * 4,
        })
        values = np.array([
            [8.0, 8.4, 7.8, 8.2, 6.0, 6.1, 6.2, 6.0],
            [4.0, 4.1, 3.9, 4.2, 5.0, 5.2, 5.1, 4.9],
            [6.0, 6.3, 5.8, 6.1, 6.0, 6.1, 5.9, 6.2],
        ])
        self.expression = pd.DataFrame(values, index=["p1", "p2", "p3"],
                                       columns=self.metadata["sample_accession"])

    def test_pinned_source_sizes(self):
        self.assertEqual((MATRIX_BYTES, ANNOTATION_BYTES), (21_530_998, 8_471_521))

    def test_pinned_md5_values(self):
        self.assertRegex(MATRIX_MD5, r"^[0-9a-f]{32}$")
        self.assertRegex(ANNOTATION_MD5, r"^[0-9a-f]{32}$")

    def test_md5_function(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "value.txt"
            path.write_bytes(b"genomics")
            self.assertEqual(md5sum(path), "e7d5b92c4f0854f9bd64ed8f6384869f")

    def test_bh_known_values(self):
        observed = benjamini_hochberg(np.array([0.01, 0.04, 0.03, 0.002]))
        np.testing.assert_allclose(observed, [0.02, 0.04, 0.04, 0.008])

    def test_bh_rejects_invalid_values(self):
        with self.assertRaises(ValueError):
            benjamini_hochberg(np.array([0.1, np.nan]))

    def test_pair_indices_match_participants(self):
        tumor, normal, patients = pair_indices(self.metadata)
        self.assertEqual((tumor, normal), ([0, 1, 2, 3], [4, 5, 6, 7]))
        self.assertEqual(patients, ["1", "2", "3", "4"])

    def test_excluded_participant_is_removed(self):
        _, _, patients = pair_indices(self.metadata, {"2"})
        self.assertEqual(patients, ["1", "3", "4"])

    def test_paired_test_returns_every_probe(self):
        result, patients = paired_test(self.expression, self.metadata)
        self.assertEqual(len(result), 3)
        self.assertEqual(len(patients), 4)
        self.assertTrue(np.isfinite(result[["p_value", "fdr"]]).all().all())

    def test_paired_effect_direction(self):
        result, _ = paired_test(self.expression, self.metadata)
        self.assertGreater(result.loc["p1", "mean_log2_fold_change"], 0)
        self.assertLess(result.loc["p2", "mean_log2_fold_change"], 0)

    def test_unambiguous_gene_symbol_rule(self):
        self.assertTrue(unambiguous_symbol("AGER"))
        self.assertFalse(unambiguous_symbol("MIR4640///DDR1"))
        self.assertFalse(unambiguous_symbol(""))

    def test_full_input_validator_rejects_wrong_dimensions(self):
        with self.assertRaisesRegex(ValueError, "Unexpected expression dimensions"):
            validate_inputs(self.expression, self.metadata)


if __name__ == "__main__":
    unittest.main(verbosity=2)
