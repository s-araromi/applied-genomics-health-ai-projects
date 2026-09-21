import io
import sys
import unittest
from pathlib import Path

PROJECT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT))

from analyze_fastq_quality import (  # noqa: E402
    FastqError,
    analyze_fastq,
    combine_metrics,
    iter_fastq_records,
    normalized_read_id,
    qc_status,
    validate_pairs,
)


FIXTURES = Path(__file__).parent / "fixtures"


class FastqQualityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.r1 = analyze_fastq(FIXTURES / "valid_R1.fastq")
        cls.r2 = analyze_fastq(FIXTURES / "valid_R2.fastq")

    def test_valid_fastq_has_two_reads(self):
        self.assertEqual(self.r1["reads"], 2)

    def test_base_and_n_counts_are_correct(self):
        self.assertEqual(self.r1["total_bases"], 10)
        self.assertEqual(self.r1["n_bases"], 1)
        self.assertAlmostEqual(self.r1["gc_percent"], 6 / 9 * 100)

    def test_phred33_quality_is_decoded(self):
        self.assertAlmostEqual(self.r1["mean_quality"], 30.0)
        self.assertAlmostEqual(self.r1["q30_percent"], 70.0)

    def test_pair_suffixes_are_normalized(self):
        self.assertEqual(normalized_read_id("PAIR001/1"), "PAIR001")
        self.assertEqual(normalized_read_id("PAIR001/2 extra"), "PAIR001")

    def test_fixture_pairs_match(self):
        validate_pairs(self.r1, self.r2)

    def test_mismatched_pair_is_rejected(self):
        altered = dict(self.r2)
        altered["identifiers"] = ["different", "PAIR002"]
        with self.assertRaises(FastqError):
            validate_pairs(self.r1, altered)

    def test_incomplete_record_is_rejected(self):
        with self.assertRaises(FastqError):
            list(iter_fastq_records(io.StringIO("@x\nACGT\n+\n")))

    def test_sequence_quality_length_mismatch_is_rejected(self):
        with self.assertRaises(FastqError):
            list(iter_fastq_records(io.StringIO("@x\nACGT\n+\nIII\n")))

    def test_invalid_base_is_rejected(self):
        with self.assertRaises(FastqError):
            list(iter_fastq_records(io.StringIO("@x\nACGX\n+\nIIII\n")))

    def test_combined_metrics_and_review_flag(self):
        combined = combine_metrics(self.r1, self.r2)
        self.assertEqual(combined["reads"], 4)
        self.assertEqual(combined["total_bases"], 20)
        self.assertEqual(qc_status(combined), "REVIEW")


if __name__ == "__main__":
    unittest.main(verbosity=2)
