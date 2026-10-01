import csv
import json
import re
import sys
import unittest
from datetime import date
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT))

from clinvar_conflicts import (  # noqa: E402
    DEFAULT_QUERY,
    SNAPSHOT_DATE,
    build_classification_counts,
    parse_conflict_explanation,
    parse_vcv_record,
    summarize_records,
)


class ClinVarConflictTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        fixture = PROJECT / "tests" / "fixtures" / "VCV000035835.90.reduced.xml"
        cls.xml = fixture.read_bytes()
        cls.record = parse_vcv_record(cls.xml, SNAPSHOT_DATE)
        with (PROJECT / "outputs" / "cftr_conflicting_variants.csv").open(
            newline="", encoding="utf-8"
        ) as handle:
            cls.output_records = list(csv.DictReader(handle))
        cls.manifest = json.loads((PROJECT / "data" / "source_manifest.json").read_text())

    def test_pinned_query_names_gene_and_conflict(self):
        self.assertEqual(
            DEFAULT_QUERY,
            'CFTR[gene] AND "conflicting classifications of pathogenicity"[clinsig]',
        )

    def test_authentic_fixture_identity(self):
        self.assertEqual(self.record["variation_id"], "35835")
        self.assertEqual(self.record["accession_version"], "VCV000035835.90")
        self.assertEqual(self.record["record_status"], "current")

    def test_authentic_fixture_grch38_fields(self):
        self.assertEqual(self.record["grch38_accession"], "NC_000007.14")
        self.assertEqual(self.record["grch38_position"], "117592169")
        self.assertEqual((self.record["grch38_ref"], self.record["grch38_alt"]), ("C", "T"))

    def test_authentic_fixture_conflict_counts(self):
        self.assertEqual(self.record["submission_count"], 31)
        self.assertEqual(self.record["pathogenic_tier_count"], 2)
        self.assertEqual(self.record["vus_count"], 16)
        self.assertEqual(self.record["benign_tier_count"], 7)
        self.assertTrue(self.record["polarized_pathogenic_benign"])

    def test_authentic_fixture_priority(self):
        self.assertEqual(self.record["priority_tier"], "Tier 1: polarized, multi-submitter")
        self.assertTrue(self.record["criteria_provided"])

    def test_conflict_explanation_parser(self):
        observed = parse_conflict_explanation(
            "Pathogenic (2); Uncertain significance (16); Benign (3); Likely benign (4)"
        )
        self.assertEqual(
            observed,
            {"Pathogenic": 2, "Uncertain significance": 16, "Benign": 3, "Likely benign": 4},
        )

    def test_malformed_conflict_explanation_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "Unexpected ClinVar conflict explanation"):
            parse_conflict_explanation("Pathogenic: two")

    def test_previous_record_status_is_rejected(self):
        altered = self.xml.replace(b"<RecordStatus>current</RecordStatus>", b"<RecordStatus>previous</RecordStatus>", 1)
        with self.assertRaisesRegex(ValueError, "not current"):
            parse_vcv_record(altered, date(2026, 9, 28))

    def test_manifest_matches_committed_record_versions(self):
        observed = [record["accession_version"] for record in self.output_records]
        self.assertEqual(observed, self.manifest["selected_accession_versions"])
        self.assertEqual(len(observed), 25)

    def test_all_committed_accessions_are_versioned_and_current(self):
        for record in self.output_records:
            self.assertRegex(record["accession_version"], r"^VCV\d{9}\.\d+$")
            self.assertEqual(record["record_status"], "current")
            self.assertEqual(record["grch38_accession"], "NC_000007.14")

    def test_classification_groups_sum_to_submission_count(self):
        for record in self.output_records:
            grouped = sum(
                int(record[column])
                for column in (
                    "pathogenic_tier_count",
                    "vus_count",
                    "benign_tier_count",
                    "other_classification_count",
                )
            )
            self.assertLessEqual(grouped, int(record["submission_count"]))

    def test_long_table_has_four_rows_per_record(self):
        rows = build_classification_counts([self.record])
        self.assertEqual(len(rows), 4)
        self.assertEqual(sum(row["submission_count"] for row in rows), 25)

    def test_summary_of_authentic_fixture(self):
        summary = summarize_records([self.record])
        self.assertEqual(summary["current_records"], 1)
        self.assertEqual(summary["polarized_records"], 1)
        self.assertEqual(summary["maximum_submissions"], 31)


if __name__ == "__main__":
    unittest.main(verbosity=2)
