"""Verify the exercise against a documented real NCBI MTHFR sequence subset."""

import unittest

from analyze_gc_content import analyze_sequence, classify_gc_content, parse_fasta


# These are genuine bases 1-120 of human MTHFR RefSeq NM_005957.5, retrieved
# from the following official, reproducible NCBI EFetch URL:
# https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NM_005957.5&seq_start=1&seq_stop=120&rettype=fasta&retmode=text
REAL_MTHFR_FIRST_120_BASES = (
    "ACTGCGTTCCCCGCCCCTGCAGCGGCCACAGTGGTGCGGCCGGCGGCCGAGCGTTCTGAGTCACCCGGGA"
    "CTGGAGGTAGGAACCCAGCCATGGTGAACGAAGCCAGAGGAAACAGCAGC"
)


class GCContentTests(unittest.TestCase):
    """Check analysis behavior without needing another network request."""

    def test_real_mthfr_subset_has_documented_gc_content(self):
        result = analyze_sequence(
            "MTHFR", "NM_005957.5", REAL_MTHFR_FIRST_120_BASES
        )

        self.assertEqual(result["length_nt"], 120)
        self.assertEqual(result["a_count"], 25)
        self.assertEqual(result["c_count"], 39)
        self.assertEqual(result["g_count"], 42)
        self.assertEqual(result["t_count"], 14)
        self.assertEqual(result["ambiguous_base_count"], 0)
        self.assertEqual(result["gc_percent"], 67.5)
        self.assertEqual(result["gc_category"], "High GC")

    def test_parser_reconstructs_real_multiline_mthfr_sequence(self):
        real_fasta = (
            ">NM_005957.5 Homo sapiens MTHFR real reference subset\n"
            f"{REAL_MTHFR_FIRST_120_BASES[:70]}\n"
            f"{REAL_MTHFR_FIRST_120_BASES[70:]}\n"
        )

        records = parse_fasta(real_fasta)

        self.assertEqual(
            records["NM_005957.5"]["sequence"], REAL_MTHFR_FIRST_120_BASES
        )

    def test_descriptive_gc_boundaries_are_consistent(self):
        self.assertEqual(classify_gc_content(39.99), "Low GC")
        self.assertEqual(classify_gc_content(40), "Moderate GC")
        self.assertEqual(classify_gc_content(60), "Moderate GC")
        self.assertEqual(classify_gc_content(60.01), "High GC")

    def test_empty_sequence_is_rejected(self):
        with self.assertRaises(ValueError):
            analyze_sequence("MTHFR", "NM_005957.5", "")


if __name__ == "__main__":
    unittest.main()
