#!/usr/bin/env python3
"""Audit high-submission CFTR records with conflicting ClinVar classifications."""

from __future__ import annotations

import argparse
import csv
import json
import time
from collections import Counter
from datetime import date
from pathlib import Path

import matplotlib.pyplot as plt

from clinvar_conflicts import (
    CLINVAR_RECORD_BASE,
    DEFAULT_QUERY,
    SNAPSHOT_DATE,
    TOP_N,
    build_classification_counts,
    fetch_vcv_xml,
    parse_vcv_record,
    retrieve_selected_subset,
    summarize_records,
    write_csv,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--email", help="Email sent to NCBI E-utilities as a contact parameter")
    parser.add_argument("--top-n", type=int, default=TOP_N, help="Number of high-submission records to audit")
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"))
    parser.add_argument("--data-dir", type=Path, default=Path("data"))
    parser.add_argument(
        "--snapshot-date",
        type=date.fromisoformat,
        default=SNAPSHOT_DATE,
        help="Snapshot date in YYYY-MM-DD format",
    )
    return parser.parse_args()


def retrieve_and_parse(
    email: str | None, top_n: int, snapshot_date: date
) -> tuple[int, int, list[dict]]:
    selection = retrieve_selected_subset(email=email, top_n=top_n)
    records = []
    for index, summary in enumerate(selection.selected_summaries, start=1):
        variation_id = str(summary["uid"])
        print(f"Validating {index}/{len(selection.selected_summaries)}: Variation ID {variation_id}")
        record = parse_vcv_record(fetch_vcv_xml(variation_id, email), snapshot_date)
        if record["accession_version"] != summary["accession_version"]:
            raise ValueError(
                f"Version mismatch for Variation ID {variation_id}: "
                f"{summary['accession_version']} versus {record['accession_version']}"
            )
        records.append(record)
        time.sleep(0.35)
    records.sort(key=lambda item: (-item["submission_count"], item["accession_version"]))
    return selection.query_count, selection.aggregate_conflict_count, records


def write_manifest(
    path: Path,
    query_count: int,
    aggregate_conflict_count: int,
    records: list[dict],
    snapshot_date: date,
) -> None:
    manifest = {
        "project": "Project 07: CFTR ClinVar classification-conflict audit",
        "snapshot_date": snapshot_date.isoformat(),
        "database": "NCBI ClinVar",
        "query": DEFAULT_QUERY,
        "query_result_count": query_count,
        "records_with_aggregate_conflicting_classification": aggregate_conflict_count,
        "subset_rule": (
            "Among query results whose ESummary aggregate germline classification was exactly "
            f"'Conflicting classifications of pathogenicity', select the top {len(records)} current records "
            "by number of supporting SCV accessions; "
            "ties resolved by versioned VCV accession"
        ),
        "selected_record_count": len(records),
        "all_selected_records_verified_current": all(record["record_status"] == "current" for record in records),
        "selected_accession_versions": [record["accession_version"] for record in records],
        "source_pages": [record["clinvar_url"] for record in records],
        "eutils_documentation": "https://www.ncbi.nlm.nih.gov/clinvar/docs/programmatic_access/",
        "review_status_documentation": "https://www.ncbi.nlm.nih.gov/clinvar/docs/review_status/",
        "usage_policy": "https://www.ncbi.nlm.nih.gov/home/about/policies/",
        "reproducibility_note": (
            "ClinVar is updated regularly. The committed CSV files preserve the parsed 2026-09-28 "
            "snapshot; a live rerun may return newer versions or a different top-25 subset."
        ),
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")


def write_summary_table(path: Path, summary: dict) -> None:
    rows = [
        {"metric": "Selected current VCV records", "value": summary["current_records"]},
        {"metric": "Records with current GRCh38 coordinates", "value": summary["records_with_current_grch38_location"]},
        {"metric": "Polarized pathogenic-versus-benign conflicts", "value": summary["polarized_records"]},
        {"metric": "Records with criteria provided", "value": summary["criteria_provided_records"]},
        {"metric": "Records last evaluated more than five years ago", "value": summary["stale_records"]},
        {"metric": "Median aggregate submissions", "value": summary["median_submissions"]},
        {"metric": "Maximum aggregate submissions", "value": summary["maximum_submissions"]},
    ]
    write_csv(path, rows)


def create_figure(records: list[dict], output_dir: Path) -> None:
    selected = records[:15]
    labels = [record["accession_version"].split(".")[0].replace("VCV000", "VCV") for record in selected]
    pathogenic = [record["pathogenic_tier_count"] for record in selected]
    uncertain = [record["vus_count"] for record in selected]
    benign = [record["benign_tier_count"] for record in selected]
    other = [record["other_classification_count"] for record in selected]

    plt.style.use("seaborn-v0_8-whitegrid")
    figure, axes = plt.subplots(1, 2, figsize=(14, 7), gridspec_kw={"width_ratios": [1.55, 1]})
    y_positions = list(range(len(selected)))
    axes[0].barh(y_positions, pathogenic, color="#C84C4C", label="Pathogenic tier")
    axes[0].barh(y_positions, uncertain, left=pathogenic, color="#D7A62B", label="Uncertain significance")
    left_benign = [a + b for a, b in zip(pathogenic, uncertain)]
    axes[0].barh(y_positions, benign, left=left_benign, color="#3B8C75", label="Benign tier")
    left_other = [a + b for a, b in zip(left_benign, benign)]
    axes[0].barh(y_positions, other, left=left_other, color="#8A8F98", label="Other")
    axes[0].set_yticks(y_positions, labels)
    axes[0].invert_yaxis()
    axes[0].set_xlabel("Submitted classifications represented in aggregate conflict")
    axes[0].set_title("Highest-submission conflicting CFTR records")
    axes[0].legend(frameon=False, fontsize=8, loc="lower right")

    tier_counts = Counter(record["priority_tier"].split(":")[0] for record in records)
    tier_order = ["Tier 1", "Tier 2", "Tier 3"]
    tier_values = [tier_counts.get(tier, 0) for tier in tier_order]
    axes[1].bar(tier_order, tier_values, color=["#355C7D", "#6C8EAD", "#B8C5D1"])
    axes[1].set_ylabel("Records")
    axes[1].set_title("Manual-curation priority tiers")
    axes[1].set_ylim(0, max(tier_values + [1]) * 1.2)
    for index, value in enumerate(tier_values):
        axes[1].text(index, value + 0.2, str(value), ha="center", va="bottom")

    figure.suptitle("Project 07 — ClinVar CFTR classification-conflict audit", fontsize=15, fontweight="bold")
    figure.text(
        0.5,
        0.01,
        "Priority tiers support evidence review only; they do not determine clinical significance.",
        ha="center",
        fontsize=9,
    )
    figure.tight_layout(rect=(0, 0.04, 1, 0.95))
    figure.savefig(output_dir / "clinvar_conflict_audit.png", dpi=220, bbox_inches="tight")
    figure.savefig(output_dir / "clinvar_conflict_audit.svg", bbox_inches="tight")
    plt.close(figure)


def write_report(
    path: Path,
    query_count: int,
    aggregate_conflict_count: int,
    records: list[dict],
    summary: dict,
    snapshot_date: date,
) -> None:
    top = records[:5]
    tier_lines = [f"- {tier}: **{count}**" for tier, count in summary["tier_counts"].items()]
    top_lines = [
        f"- [{record['accession_version']}]({record['clinvar_url']}): "
        f"{record['submission_count']} submissions from {record['submitter_count']} submitters; "
        f"{record['priority_tier']}."
        for record in top
    ]
    text = [
        "# CFTR ClinVar Classification-Conflict Audit",
        "",
        f"**Snapshot date:** {snapshot_date.isoformat()}",
        "",
        "## Scope",
        "",
        f"The exact ClinVar query returned **{query_count}** records; **{aggregate_conflict_count}** had the exact "
        "aggregate germline label `Conflicting classifications of pathogenicity` in ESummary. The audit examined "
        f"the **{len(records)}** eligible records with the largest numbers of supporting SCV accessions. Every selected VCV archive resolved "
        "as `current` during retrieval.",
        "",
        "## Verified results",
        "",
        f"- Current records: **{summary['current_records']}/{len(records)}**",
        f"- Records with a current GRCh38 coordinate: **{summary['records_with_current_grch38_location']}/{len(records)}**",
        f"- Polarized pathogenic-tier versus benign-tier disagreements: **{summary['polarized_records']}**",
        f"- Records with criteria provided: **{summary['criteria_provided_records']}**",
        f"- Records last evaluated more than five years before the snapshot: **{summary['stale_records']}**",
        f"- Median aggregate classification count: **{summary['median_submissions']:g}**",
        f"- Maximum aggregate classification count: **{summary['maximum_submissions']}**",
        "",
        "## Priority tiers",
        "",
        *tier_lines,
        "",
        "Tier 1 requires a pathogenic-tier/benign-tier polarization, criteria-provided review status, and at "
        "least three submitters. Tier 2 requires criteria-provided review status and at least three submitters "
        "without that polarization. Tier 3 contains records with more limited review support. These are "
        "workflow priorities for manual evidence reconciliation, not variant classifications.",
        "",
        "## Highest-submission records",
        "",
        *top_lines,
        "",
        "## Interpretation",
        "",
        "A large submission count does not make one interpretation correct. It indicates a larger evidence-"
        "reconciliation workload. Polarized records deserve particular attention because the aggregate conflict "
        "crosses the pathogenic-to-benign boundary. Clinical interpretation still requires condition-specific "
        "review of individual SCV assertions, inheritance, phenotype, allele frequency, functional evidence, and "
        "the most recent expert guidance.",
        "",
        "## Limitations",
        "",
        "- This is a high-submission subset, not all conflicting CFTR records.",
        "- ClinVar is dynamic; later reruns can change record versions and rankings.",
        "- Aggregate counts can combine condition-specific assertions and should not be treated as votes.",
        "- The priority tiers are transparent operational rules, not validated clinical decision rules.",
        "- No patient-level data were used, and the outputs must not be used for diagnosis or treatment.",
        "",
        "## Source",
        "",
        f"- [ClinVar query](https://www.ncbi.nlm.nih.gov/clinvar/?term={DEFAULT_QUERY.replace(' ', '+')})",
        "- [ClinVar programmatic access](https://www.ncbi.nlm.nih.gov/clinvar/docs/programmatic_access/)",
        "- [ClinVar review status](https://www.ncbi.nlm.nih.gov/clinvar/docs/review_status/)",
        "",
    ]
    path.write_text("\n".join(text), encoding="utf-8")


def main() -> None:
    args = parse_args()
    if args.top_n < 1 or args.top_n > 100:
        raise ValueError("--top-n must be between 1 and 100")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.data_dir.mkdir(parents=True, exist_ok=True)
    query_count, aggregate_conflict_count, records = retrieve_and_parse(
        args.email, args.top_n, args.snapshot_date
    )
    summary = summarize_records(records)

    write_csv(args.output_dir / "cftr_conflicting_variants.csv", records)
    write_csv(args.output_dir / "conflict_classification_counts.csv", build_classification_counts(records))
    write_summary_table(args.output_dir / "conflict_audit_summary.csv", summary)
    write_manifest(
        args.data_dir / "source_manifest.json",
        query_count,
        aggregate_conflict_count,
        records,
        args.snapshot_date,
    )
    create_figure(records, args.output_dir)
    write_report(
        args.output_dir / "clinvar_conflict_audit_report.md",
        query_count,
        aggregate_conflict_count,
        records,
        summary,
        args.snapshot_date,
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
