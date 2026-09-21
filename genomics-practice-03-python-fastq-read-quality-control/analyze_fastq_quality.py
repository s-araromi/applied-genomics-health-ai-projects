#!/usr/bin/env python3
"""Quality-control analysis of a reproducible paired FASTQ subset from ERR194147."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import sys
import urllib.request
from collections import Counter
from contextlib import contextmanager
from pathlib import Path

RUN_ACCESSION = "ERR194147"
SUBSET_READS = 5_000
R1_URL = "https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR194/ERR194147/ERR194147_1.fastq.gz"
R2_URL = "https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR194/ERR194147/ERR194147_2.fastq.gz"
VALID_BASES = set("ACGTN")


class FastqError(ValueError):
    """Raised when FASTQ content is structurally invalid."""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reads", type=int, default=SUBSET_READS)
    parser.add_argument("--r1", type=Path, help="Existing R1 FASTQ or FASTQ.GZ file")
    parser.add_argument("--r2", type=Path, help="Existing R2 FASTQ or FASTQ.GZ file")
    parser.add_argument("--data-dir", type=Path, default=Path("data"))
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"))
    return parser.parse_args()


@contextmanager
def open_fastq(path: Path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="ascii", newline="") as handle:
        yield handle


def iter_fastq_records(handle):
    """Yield validated (identifier, sequence, quality) FASTQ records."""
    record_number = 0
    while True:
        header = handle.readline()
        if not header:
            return
        sequence = handle.readline()
        plus = handle.readline()
        quality = handle.readline()
        record_number += 1
        if not sequence or not plus or not quality:
            raise FastqError(f"Record {record_number} is incomplete")
        header = header.rstrip("\r\n")
        sequence = sequence.rstrip("\r\n").upper()
        plus = plus.rstrip("\r\n")
        quality = quality.rstrip("\r\n")
        if not header.startswith("@"):
            raise FastqError(f"Record {record_number} header does not start with @")
        if not plus.startswith("+"):
            raise FastqError(f"Record {record_number} separator does not start with +")
        if len(sequence) != len(quality):
            raise FastqError(f"Record {record_number} sequence and quality lengths differ")
        invalid = set(sequence) - VALID_BASES
        if invalid:
            raise FastqError(f"Record {record_number} contains invalid bases: {sorted(invalid)}")
        if any(ord(char) < 33 or ord(char) > 126 for char in quality):
            raise FastqError(f"Record {record_number} has invalid quality characters")
        yield header[1:], sequence, quality


def normalized_read_id(identifier: str) -> str:
    token = identifier.split()[0]
    if token.endswith("/1") or token.endswith("/2"):
        token = token[:-2]
    return token


def download_subset(url: str, destination: Path, reads: int) -> None:
    """Stream-decompress only the leading requested records into a small gzip file."""
    if reads <= 0:
        raise ValueError("reads must be positive")
    destination.parent.mkdir(parents=True, exist_ok=True)
    request = urllib.request.Request(url, headers={"User-Agent": "genomics-practice-fastq-qc/1.0"})
    try:
        response = urllib.request.urlopen(request, timeout=120)
        source_binary = gzip.GzipFile(fileobj=response)
        source = io.TextIOWrapper(source_binary, encoding="ascii", newline="")
        output_binary = destination.open("wb")
        compressed = gzip.GzipFile(filename="", mode="wb", fileobj=output_binary, mtime=0)
        output = io.TextIOWrapper(compressed, encoding="ascii", newline="\n")
        count = 0
        try:
            for identifier, sequence, quality in iter_fastq_records(source):
                output.write(f"@{identifier}\n{sequence}\n+\n{quality}\n")
                count += 1
                if count == reads:
                    break
        finally:
            output.close()
            source.close()
            response.close()
        if count != reads:
            destination.unlink(missing_ok=True)
            raise FastqError(f"Source ended after {count:,} records; expected {reads:,}")
    except Exception:
        destination.unlink(missing_ok=True)
        raise


def analyze_fastq(path: Path) -> dict:
    reads = 0
    total_bases = 0
    q30_bases = 0
    quality_sum = 0
    base_counts = Counter()
    cycle_quality_sum: list[int] = []
    cycle_q30: list[int] = []
    cycle_depth: list[int] = []
    identifiers: list[str] = []
    lengths: list[int] = []

    with open_fastq(path) as handle:
        for identifier, sequence, quality_string in iter_fastq_records(handle):
            scores = [ord(char) - 33 for char in quality_string]
            reads += 1
            identifiers.append(normalized_read_id(identifier))
            lengths.append(len(sequence))
            total_bases += len(sequence)
            base_counts.update(sequence)
            quality_sum += sum(scores)
            q30_bases += sum(score >= 30 for score in scores)
            while len(cycle_quality_sum) < len(scores):
                cycle_quality_sum.append(0)
                cycle_q30.append(0)
                cycle_depth.append(0)
            for index, score in enumerate(scores):
                cycle_quality_sum[index] += score
                cycle_q30[index] += score >= 30
                cycle_depth[index] += 1

    if reads == 0:
        raise FastqError(f"{path} contains no FASTQ records")
    known_bases = sum(base_counts[base] for base in "ACGT")
    cycles = [
        {
            "cycle": index + 1,
            "mean_quality": quality_sum_cycle / cycle_depth[index],
            "q30_percent": 100 * cycle_q30[index] / cycle_depth[index],
        }
        for index, quality_sum_cycle in enumerate(cycle_quality_sum)
    ]
    return {
        "reads": reads,
        "total_bases": total_bases,
        "mean_quality": quality_sum / total_bases,
        "q30_percent": 100 * q30_bases / total_bases,
        "gc_percent": 100 * (base_counts["G"] + base_counts["C"]) / known_bases,
        "n_bases": base_counts["N"],
        "n_percent": 100 * base_counts["N"] / total_bases,
        "min_length": min(lengths),
        "max_length": max(lengths),
        "base_counts": base_counts,
        "identifiers": identifiers,
        "quality_sum": quality_sum,
        "q30_bases": q30_bases,
        "known_bases": known_bases,
        "cycles": cycles,
    }


def validate_pairs(r1: dict, r2: dict) -> None:
    if r1["reads"] != r2["reads"]:
        raise FastqError("R1 and R2 contain different numbers of reads")
    for index, (r1_id, r2_id) in enumerate(zip(r1["identifiers"], r2["identifiers"]), start=1):
        if r1_id != r2_id:
            raise FastqError(f"Read-pair identifiers differ at pair {index}: {r1_id} != {r2_id}")


def combine_metrics(r1: dict, r2: dict) -> dict:
    total_bases = r1["total_bases"] + r2["total_bases"]
    known_bases = r1["known_bases"] + r2["known_bases"]
    gc_bases = sum(r1["base_counts"][base] + r2["base_counts"][base] for base in "GC")
    return {
        "reads": r1["reads"] + r2["reads"],
        "total_bases": total_bases,
        "mean_quality": (r1["quality_sum"] + r2["quality_sum"]) / total_bases,
        "q30_percent": 100 * (r1["q30_bases"] + r2["q30_bases"]) / total_bases,
        "gc_percent": 100 * gc_bases / known_bases,
        "n_bases": r1["n_bases"] + r2["n_bases"],
        "n_percent": 100 * (r1["n_bases"] + r2["n_bases"]) / total_bases,
        "min_length": min(r1["min_length"], r2["min_length"]),
        "max_length": max(r1["max_length"], r2["max_length"]),
    }


def qc_status(metrics: dict) -> str:
    return "PASS" if metrics["q30_percent"] >= 90 and metrics["n_percent"] <= 0.1 else "REVIEW"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_outputs(output_dir: Path, r1_path: Path, r2_path: Path, r1: dict, r2: dict) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    combined = combine_metrics(r1, r2)
    rows = [("R1", r1), ("R2", r2), ("Combined", combined)]
    with (output_dir / "fastq_qc_summary.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["mate", "reads", "bases", "mean_phred", "q30_percent", "gc_percent", "n_bases", "n_percent", "min_length", "max_length", "status"])
        for mate, metrics in rows:
            writer.writerow([mate, metrics["reads"], metrics["total_bases"], f'{metrics["mean_quality"]:.2f}', f'{metrics["q30_percent"]:.2f}', f'{metrics["gc_percent"]:.2f}', metrics["n_bases"], f'{metrics["n_percent"]:.5f}', metrics["min_length"], metrics["max_length"], qc_status(metrics)])

    with (output_dir / "per_cycle_quality.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["mate", "cycle", "mean_phred", "q30_percent"])
        for mate, metrics in (("R1", r1), ("R2", r2)):
            for cycle in metrics["cycles"]:
                writer.writerow([mate, cycle["cycle"], f'{cycle["mean_quality"]:.4f}', f'{cycle["q30_percent"]:.4f}'])

    with (output_dir / "subset_checksums.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["file", "bytes", "sha256"])
        for path in (r1_path, r2_path):
            writer.writerow([path.name, path.stat().st_size, sha256(path)])

    report = f"""# Paired-End FASTQ Quality-Control Report

## Dataset

- ENA run: `{RUN_ACCESSION}`
- Sample: NA12878
- Subset: first {r1['reads']:,} read pairs
- Read length: {r1['min_length']}–{r1['max_length']} bases (R1); {r2['min_length']}–{r2['max_length']} bases (R2)

## Results

| Mate | Reads | Bases | Mean Phred | Q30 bases | GC content | N bases | Status |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| R1 | {r1['reads']:,} | {r1['total_bases']:,} | {r1['mean_quality']:.2f} | {r1['q30_percent']:.2f}% | {r1['gc_percent']:.2f}% | {r1['n_bases']:,} | {qc_status(r1)} |
| R2 | {r2['reads']:,} | {r2['total_bases']:,} | {r2['mean_quality']:.2f} | {r2['q30_percent']:.2f}% | {r2['gc_percent']:.2f}% | {r2['n_bases']:,} | {qc_status(r2)} |
| Combined | {combined['reads']:,} | {combined['total_bases']:,} | {combined['mean_quality']:.2f} | {combined['q30_percent']:.2f}% | {combined['gc_percent']:.2f}% | {combined['n_bases']:,} | {qc_status(combined)} |

## Interpretation

Q30 represents an estimated base-call error probability of 0.1%. A mate is marked
`REVIEW` when fewer than 90% of bases reach Q30 or more than 0.1% are ambiguous
`N` calls. A review flag is descriptive and does not by itself make the full run
unsuitable for downstream analysis.
"""
    (output_dir / "fastq_qc_report.md").write_text(report, encoding="utf-8")

    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise RuntimeError("matplotlib is required; run: python -m pip install -r requirements.txt") from exc

    figure, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    for mate, metrics, color in (("R1", r1, "#176B87"), ("R2", r2, "#D97706")):
        axes[0].plot([x["cycle"] for x in metrics["cycles"]], [x["mean_quality"] for x in metrics["cycles"]], label=mate, color=color, linewidth=2)
    axes[0].axhline(30, color="#6B7280", linestyle="--", linewidth=1)
    axes[0].set(title="Mean quality by sequencing cycle", xlabel="Cycle", ylabel="Mean Phred score")
    axes[0].legend(frameon=False)
    axes[0].grid(alpha=0.2)
    labels = ["Mean Phred", "Q30 bases (%)", "GC content (%)"]
    x = range(len(labels))
    width = 0.35
    axes[1].bar([i - width / 2 for i in x], [r1["mean_quality"], r1["q30_percent"], r1["gc_percent"]], width, label="R1", color="#176B87")
    axes[1].bar([i + width / 2 for i in x], [r2["mean_quality"], r2["q30_percent"], r2["gc_percent"]], width, label="R2", color="#D97706")
    axes[1].set(title="Read-level quality summary", xticks=list(x), xticklabels=labels)
    axes[1].legend(frameon=False)
    axes[1].grid(axis="y", alpha=0.2)
    figure.suptitle(f"{RUN_ACCESSION}: first {r1['reads']:,} paired reads", fontsize=13, fontweight="bold")
    figure.tight_layout()
    figure.savefig(output_dir / "fastq_quality_overview.png", dpi=180, bbox_inches="tight")
    figure.savefig(output_dir / "fastq_quality_overview.svg", bbox_inches="tight")
    plt.close(figure)


def main() -> int:
    args = parse_args()
    if args.reads <= 0:
        raise ValueError("--reads must be positive")
    if bool(args.r1) != bool(args.r2):
        raise ValueError("Provide both --r1 and --r2, or neither")
    if args.r1:
        r1_path, r2_path = args.r1, args.r2
    else:
        r1_path = args.data_dir / f"{RUN_ACCESSION}_R1_first{args.reads}.fastq.gz"
        r2_path = args.data_dir / f"{RUN_ACCESSION}_R2_first{args.reads}.fastq.gz"
        if not r1_path.exists():
            print(f"Streaming the first {args.reads:,} R1 records from ENA...")
            download_subset(R1_URL, r1_path, args.reads)
        if not r2_path.exists():
            print(f"Streaming the first {args.reads:,} R2 records from ENA...")
            download_subset(R2_URL, r2_path, args.reads)

    print("Validating and analyzing paired FASTQ records...")
    r1 = analyze_fastq(r1_path)
    r2 = analyze_fastq(r2_path)
    validate_pairs(r1, r2)
    write_outputs(args.output_dir, r1_path, r2_path, r1, r2)
    combined = combine_metrics(r1, r2)
    print(f"R1: {r1['reads']:,} reads | Q{r1['mean_quality']:.2f} | {r1['q30_percent']:.2f}% Q30 | {r1['gc_percent']:.2f}% GC | {r1['n_bases']:,} N")
    print(f"R2: {r2['reads']:,} reads | Q{r2['mean_quality']:.2f} | {r2['q30_percent']:.2f}% Q30 | {r2['gc_percent']:.2f}% GC | {r2['n_bases']:,} N")
    print(f"Combined: {combined['reads']:,} reads | Q{combined['mean_quality']:.2f} | {combined['q30_percent']:.2f}% Q30")
    print(f"Outputs written to {args.output_dir.resolve()}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (FastqError, OSError, RuntimeError, ValueError) as exc:
        print(f"Analysis stopped: {exc}", file=sys.stderr)
        raise SystemExit(1)
