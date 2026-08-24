#!/usr/bin/env python3
"""Analyze GC content in four real human NCBI RefSeq transcript sequences.

Run from this exercise directory:

    python analyze_gc_content.py --email your-real-email@example.com

Only Python's standard library is required.
"""

import argparse
import csv
from pathlib import Path
from statistics import mean
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen


# An accession is a database identifier. The number after the period is the
# sequence version, making the exact reference used in this project traceable.
# AGT NM_000029.4 is intentionally excluded because NCBI has suppressed it.
GENE_ACCESSIONS = {
    "MTHFR": "NM_005957.5",
    "ACE": "NM_000789.4",
    "AGT": "NM_001382817.3",
    "NOS3": "NM_000603.5",
}

# EFetch is NCBI's official service for downloading database records.
EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

# Resolve paths relative to this script, so output files go into the correct
# exercise folder even if the command is started from another directory.
PROJECT_DIRECTORY = Path(__file__).resolve().parent
OUTPUT_DIRECTORY = PROJECT_DIRECTORY / "outputs"


def build_download_url(email=None):
    """Construct one NCBI request for all four real RefSeq records."""
    parameters = {
        "db": "nuccore",                  # NCBI's nucleotide database.
        "id": ",".join(GENE_ACCESSIONS.values()),
        "rettype": "fasta",              # Return sequences in FASTA format.
        "retmode": "text",              # Return ordinary readable text.
        "tool": "applied_genomics_health_ai_projects",
    }

    # NCBI recommends including a valid contact address with automated requests.
    # It is supplied at runtime instead of being committed to a public repository.
    if email:
        parameters["email"] = email

    # URL encoding safely converts a dictionary into an HTTP query string.
    return f"{EFETCH_URL}?{urlencode(parameters)}"


def download_fasta(email=None):
    """Download genuine sequence records directly from NCBI."""
    request = Request(
        build_download_url(email),
        headers={"User-Agent": "applied-genomics-health-ai-projects/1.0"},
    )

    # The timeout prevents the program from waiting forever if NCBI is slow.
    with urlopen(request, timeout=60) as response:
        fasta_text = response.read().decode("utf-8")

    if not fasta_text.lstrip().startswith(">"):
        raise ValueError("NCBI did not return recognizable FASTA data.")

    return fasta_text


def parse_fasta(fasta_text):
    """Convert multi-record FASTA text into accession-to-sequence entries."""
    records = {}
    current_accession = None
    current_description = ""
    sequence_lines = []

    for raw_line in fasta_text.splitlines():
        line = raw_line.strip()

        # Ignore blank lines between records.
        if not line:
            continue

        if line.startswith(">"):
            # Before starting a new record, save the previous record, if any.
            if current_accession is not None:
                records[current_accession] = {
                    "description": current_description,
                    "sequence": "".join(sequence_lines).upper(),
                }

            # A FASTA header starts with ">"; its first word is the accession.
            header = line[1:]
            current_accession = header.split(maxsplit=1)[0]
            current_description = header
            sequence_lines = []
        else:
            if current_accession is None:
                raise ValueError("A sequence appeared before its FASTA header.")
            sequence_lines.append(line)

    # The final record has no subsequent header, so save it after the loop ends.
    if current_accession is not None:
        records[current_accession] = {
            "description": current_description,
            "sequence": "".join(sequence_lines).upper(),
        }

    if not records:
        raise ValueError("No FASTA records could be parsed.")

    return records


def classify_gc_content(gc_percent):
    """Apply illustrative, non-clinical sequence-composition categories."""
    if gc_percent < 40:
        return "Low GC"
    if gc_percent <= 60:
        return "Moderate GC"
    return "High GC"


def analyze_sequence(gene, accession, sequence):
    """Count bases and calculate GC percentage for one genuine transcript."""
    sequence = sequence.upper()

    # NCBI transcript FASTA uses the DNA-style alphabet A/C/G/T even though the
    # record describes mRNA; this is a normal database representation.
    base_counts = {base: sequence.count(base) for base in "ACGT"}
    known_base_count = sum(base_counts.values())
    ambiguous_base_count = len(sequence) - known_base_count

    if known_base_count == 0:
        raise ValueError(f"{gene} has no recognized A, C, G, or T bases.")

    gc_count = base_counts["G"] + base_counts["C"]
    gc_percent = round((gc_count / known_base_count) * 100, 2)

    return {
        "gene": gene,
        "accession": accession,
        "length_nt": len(sequence),
        "a_count": base_counts["A"],
        "c_count": base_counts["C"],
        "g_count": base_counts["G"],
        "t_count": base_counts["T"],
        "ambiguous_base_count": ambiguous_base_count,
        "gc_count": gc_count,
        "gc_percent": gc_percent,
        "gc_category": classify_gc_content(gc_percent),
    }


def analyze_records(records):
    """Check that NCBI returned every requested accession, then analyze it."""
    missing_accessions = [
        accession
        for accession in GENE_ACCESSIONS.values()
        if accession not in records
    ]

    if missing_accessions:
        missing_text = ", ".join(missing_accessions)
        raise ValueError(f"NCBI did not return these requested records: {missing_text}")

    results = []

    for gene, accession in GENE_ACCESSIONS.items():
        sequence = records[accession]["sequence"]
        results.append(analyze_sequence(gene, accession, sequence))

    return results


def write_csv(results):
    """Save a spreadsheet-compatible table of all analysis results."""
    csv_path = OUTPUT_DIRECTORY / "gc_content_results.csv"

    with csv_path.open("w", newline="", encoding="utf-8") as output_file:
        writer = csv.DictWriter(output_file, fieldnames=results[0].keys())
        writer.writeheader()
        writer.writerows(results)

    return csv_path


def write_markdown_report(results):
    """Save a human-readable summary suitable for display on GitHub."""
    report_path = OUTPUT_DIRECTORY / "gc_content_report.md"
    ranked_results = sorted(results, key=lambda row: row["gc_percent"], reverse=True)
    highest = ranked_results[0]
    lowest = ranked_results[-1]

    # Each transcript contributes equally to the unweighted mean.
    average_gc = mean(row["gc_percent"] for row in results)

    # The pooled value combines all counted bases; longer transcripts contribute
    # more strongly, which can produce a different overall GC percentage.
    total_gc_bases = sum(row["gc_count"] for row in results)
    total_known_bases = sum(
        row["length_nt"] - row["ambiguous_base_count"]
        for row in results
    )
    pooled_gc = (total_gc_bases / total_known_bases) * 100

    lines = [
        "# GC-content results: real NCBI RefSeq transcripts",
        "",
        "| Rank | Gene | RefSeq accession | Length (nt) | GC (%) | Category |",
        "| --- | --- | --- | ---: | ---: | --- |",
    ]

    for rank, row in enumerate(ranked_results, start=1):
        lines.append(
            f"| {rank} | {row['gene']} | {row['accession']} | "
            f"{row['length_nt']:,} | {row['gc_percent']:.2f} | "
            f"{row['gc_category']} |"
        )

    lines.extend([
        "",
        f"- Highest GC content: **{highest['gene']} ({highest['gc_percent']:.2f}%)**.",
        f"- Lowest GC content: **{lowest['gene']} ({lowest['gc_percent']:.2f}%)**.",
        f"- Mean transcript-level GC content: **{average_gc:.2f}%**.",
        f"- Pooled, length-weighted GC content: **{pooled_gc:.2f}%**.",
        "",
        "## Interpretation caution",
        "",
        "These values describe full reference-transcript composition only. They do",
        "not measure gene expression, establish disease risk, diagnose a patient,",
        "or determine PCR primer performance. The 40% and 60% labels are",
        "illustrative educational categories, not clinical cutoffs.",
        "",
        "## Data source",
        "",
        "NCBI Nucleotide / RefSeq: https://www.ncbi.nlm.nih.gov/refseq/",
        "",
    ])

    report_path.write_text("\n".join(lines), encoding="utf-8")
    return report_path


def print_console_summary(results):
    """Display the most important findings in the terminal."""
    print("\nGC-content results for real NCBI RefSeq transcripts\n")

    for row in results:
        print(
            f"{row['gene']:5} | {row['accession']:14} | "
            f"{row['length_nt']:5,} nt | "
            f"GC: {row['gc_percent']:6.2f}% | {row['gc_category']}"
        )

    highest = max(results, key=lambda row: row["gc_percent"])
    lowest = min(results, key=lambda row: row["gc_percent"])
    average_gc = mean(row["gc_percent"] for row in results)

    print(f"\nHighest GC content: {highest['gene']} ({highest['gc_percent']:.2f}%)")
    print(f"Lowest GC content: {lowest['gene']} ({lowest['gc_percent']:.2f}%)")
    print(f"Mean transcript-level GC content: {average_gc:.2f}%")


def main():
    """Run the complete download, quality-control, analysis, and export workflow."""
    parser = argparse.ArgumentParser(
        description="Analyze GC content in four real human NCBI RefSeq transcripts."
    )
    parser.add_argument(
        "--email",
        help="Your valid email address; NCBI recommends including it with API requests.",
    )
    arguments = parser.parse_args()

    if not arguments.email:
        print("Note: NCBI recommends adding --email your-real-email@example.com")

    print("Downloading four real, versioned RefSeq transcripts from NCBI...")

    try:
        fasta_text = download_fasta(arguments.email)
        records = parse_fasta(fasta_text)
        results = analyze_records(records)
    except HTTPError as error:
        raise SystemExit(f"NCBI returned HTTP error {error.code}: {error.reason}")
    except URLError as error:
        raise SystemExit(f"Could not connect to NCBI: {error.reason}")
    except (TimeoutError, ValueError) as error:
        raise SystemExit(f"Analysis stopped: {error}")

    OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)
    fasta_path = OUTPUT_DIRECTORY / "refseq_transcripts.fasta"
    fasta_path.write_text(fasta_text, encoding="utf-8")

    csv_path = write_csv(results)
    report_path = write_markdown_report(results)
    print_console_summary(results)

    print("\nFiles created:")
    print(f"- {fasta_path}")
    print(f"- {csv_path}")
    print(f"- {report_path}")


# This guard runs main() only when the file is executed directly. It prevents a
# download when the test suite imports individual functions for verification.
if __name__ == "__main__":
    main()
