"""Retrieval, parsing, validation, and reporting functions for Project 07."""

from __future__ import annotations

import csv
import json
import math
import re
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from collections import Counter
from dataclasses import dataclass
from datetime import date, datetime
from pathlib import Path
from typing import Any, Iterable


EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
CLINVAR_RECORD_BASE = "https://www.ncbi.nlm.nih.gov/clinvar/variation"
DEFAULT_QUERY = 'CFTR[gene] AND "conflicting classifications of pathogenicity"[clinsig]'
GENE_SYMBOL = "CFTR"
TOP_N = 25
SNAPSHOT_DATE = date(2026, 9, 28)
USER_AGENT = "applied-genomics-health-ai-projects/07-clinvar-conflict-audit"


@dataclass(frozen=True)
class RetrievalResult:
    """Selected ClinVar summaries and provenance for a live query."""

    query_count: int
    aggregate_conflict_count: int
    selected_summaries: list[dict[str, Any]]
    selected_ids: list[str]


def _request_bytes(
    endpoint: str,
    parameters: dict[str, Any],
    *,
    email: str | None = None,
    attempts: int = 4,
    timeout: int = 90,
) -> bytes:
    """Call one NCBI E-utilities endpoint with bounded retries."""
    payload = dict(parameters)
    payload["tool"] = "clinvar_conflict_audit"
    if email:
        payload["email"] = email
    url = f"{EUTILS_BASE}/{endpoint}?{urllib.parse.urlencode(payload)}"
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    for attempt in range(attempts):
        try:
            with urllib.request.urlopen(request, timeout=timeout) as response:
                return response.read()
        except (urllib.error.URLError, TimeoutError) as error:
            if attempt == attempts - 1:
                raise RuntimeError(f"NCBI request failed after {attempts} attempts: {url}") from error
            time.sleep(2**attempt)
    raise AssertionError("unreachable")


def _request_json(endpoint: str, parameters: dict[str, Any], email: str | None = None) -> dict[str, Any]:
    return json.loads(_request_bytes(endpoint, parameters, email=email))


def search_clinvar(query: str = DEFAULT_QUERY, email: str | None = None) -> tuple[int, list[str]]:
    """Return the complete variation-ID list for a ClinVar query."""
    initial = _request_json(
        "esearch.fcgi",
        {"db": "clinvar", "term": query, "retmode": "json", "retmax": 0},
        email,
    )["esearchresult"]
    count = int(initial["count"])
    if count == 0:
        raise ValueError("The ClinVar query returned no records")
    complete = _request_json(
        "esearch.fcgi",
        {"db": "clinvar", "term": query, "retmode": "json", "retmax": count},
        email,
    )["esearchresult"]
    ids = complete["idlist"]
    if len(ids) != count or len(ids) != len(set(ids)):
        raise ValueError("ClinVar search identifiers were incomplete or duplicated")
    return count, ids


def retrieve_summaries(ids: list[str], email: str | None = None, batch_size: int = 100) -> list[dict[str, Any]]:
    """Retrieve ESummary documents in URL-safe batches."""
    summaries: list[dict[str, Any]] = []
    for start in range(0, len(ids), batch_size):
        batch = ids[start : start + batch_size]
        response = _request_json(
            "esummary.fcgi",
            {"db": "clinvar", "id": ",".join(batch), "retmode": "json"},
            email,
        )["result"]
        missing = [uid for uid in batch if uid not in response]
        if missing:
            raise ValueError(f"Missing ClinVar summaries: {', '.join(missing[:5])}")
        summaries.extend(response[uid] for uid in batch)
        time.sleep(0.35)
    return summaries


def _summary_scv_count(summary: dict[str, Any]) -> int:
    return len(summary.get("supporting_submissions", {}).get("scv", []))


def validate_summary(summary: dict[str, Any]) -> None:
    """Reject summaries that do not meet the pre-specified source checks."""
    accession = summary.get("accession_version", "")
    if not re.fullmatch(r"VCV\d{9}\.\d+", accession):
        raise ValueError(f"Missing versioned VCV accession for UID {summary.get('uid', '')}")
    genes = {entry.get("symbol") for entry in summary.get("genes", [])}
    if GENE_SYMBOL not in genes:
        raise ValueError(f"{accession} is not assigned to {GENE_SYMBOL}")
    classification = summary.get("germline_classification", {})
    if classification.get("description") != "Conflicting classifications of pathogenicity":
        raise ValueError(f"{accession} no longer has the requested aggregate classification")


def select_high_submission_subset(
    summaries: list[dict[str, Any]], top_n: int = TOP_N
) -> list[dict[str, Any]]:
    """Select a deterministic, genuine subset with the largest SCV counts."""
    eligible = [
        summary
        for summary in summaries
        if summary.get("germline_classification", {}).get("description")
        == "Conflicting classifications of pathogenicity"
    ]
    for summary in eligible:
        validate_summary(summary)
    ranked = sorted(
        eligible,
        key=lambda item: (-_summary_scv_count(item), item["accession_version"]),
    )
    return ranked[: min(top_n, len(ranked))]


def retrieve_selected_subset(
    query: str = DEFAULT_QUERY, email: str | None = None, top_n: int = TOP_N
) -> RetrievalResult:
    count, ids = search_clinvar(query, email)
    summaries = retrieve_summaries(ids, email)
    selected = select_high_submission_subset(summaries, top_n)
    aggregate_conflict_count = sum(
        summary.get("germline_classification", {}).get("description")
        == "Conflicting classifications of pathogenicity"
        for summary in summaries
    )
    return RetrievalResult(
        count, aggregate_conflict_count, selected, [str(item["uid"]) for item in selected]
    )


def fetch_vcv_xml(variation_id: str, email: str | None = None) -> bytes:
    """Fetch a current VCV archive by ClinVar Variation ID."""
    xml = _request_bytes(
        "efetch.fcgi",
        {"db": "clinvar", "id": variation_id, "rettype": "vcv", "is_variationid": "true"},
        email=email,
    )
    root = ET.fromstring(xml)
    records = root.findall("VariationArchive")
    if len(records) != 1:
        raise ValueError(f"Expected one VCV record for Variation ID {variation_id}")
    return xml


def parse_conflict_explanation(explanation: str) -> dict[str, int]:
    """Parse ClinVar's aggregate classification-count explanation."""
    counts: dict[str, int] = {}
    if not explanation.strip():
        return counts
    for component in explanation.split(";"):
        match = re.fullmatch(r"\s*(.+?)\s*\((\d+)\)\s*", component)
        if not match:
            raise ValueError(f"Unexpected ClinVar conflict explanation component: {component!r}")
        counts[match.group(1)] = int(match.group(2))
    return counts


def _unique_text(elements: Iterable[ET.Element]) -> list[str]:
    return sorted({element.text.strip() for element in elements if element.text and element.text.strip()})


def _find_grch38_location(record: ET.Element) -> ET.Element | None:
    locations = record.findall("./ClassifiedRecord/SimpleAllele/Location/SequenceLocation")
    current = [
        location
        for location in locations
        if location.get("Assembly") == "GRCh38" and location.get("AssemblyStatus") == "current"
    ]
    if not current:
        return None
    return next((location for location in current if location.get("forDisplay") == "true"), current[0])


def parse_vcv_record(xml: bytes, snapshot_date: date = SNAPSHOT_DATE) -> dict[str, Any]:
    """Parse and quality-check one current ClinVar VCV record."""
    root = ET.fromstring(xml)
    archive = root.find("VariationArchive")
    if archive is None:
        raise ValueError("ClinVar XML contains no VariationArchive")
    record_status = (archive.findtext("RecordStatus") or "").strip()
    if record_status != "current":
        raise ValueError(
            f"VCV{int(archive.get('VariationID', '0')):09d} is {record_status or 'status unknown'}, not current"
        )
    classification = archive.find("./ClassifiedRecord/Classifications/GermlineClassification")
    if classification is None:
        raise ValueError("Current ClinVar archive has no aggregate germline classification")
    description = (classification.findtext("Description") or "").strip()
    if description != "Conflicting classifications of pathogenicity":
        raise ValueError(f"Unexpected aggregate germline classification: {description}")
    explanation = (classification.findtext("Explanation") or "").strip()
    class_counts = parse_conflict_explanation(explanation)
    location = _find_grch38_location(archive)
    genes = {node.get("Symbol", "") for node in archive.findall("./ClassifiedRecord/SimpleAllele/GeneList/Gene")}
    if GENE_SYMBOL not in genes:
        raise ValueError("Current VCV record is not assigned to CFTR")

    date_text = classification.get("DateLastEvaluated", "")
    evaluated = datetime.strptime(date_text, "%Y-%m-%d").date() if date_text else None
    age_years = ((snapshot_date - evaluated).days / 365.25) if evaluated else math.nan
    traits = _unique_text(
        archive.findall("./ClassifiedRecord/RCVList/RCVAccession/ClassifiedConditionList/ClassifiedCondition")
    )
    informative_traits = [name for name in traits if name.lower() not in {"not provided", "not specified"}]
    consequences = sorted(
        {
            node.get("Type", "").strip()
            for node in archive.findall("./ClassifiedRecord/SimpleAllele/HGVSlist/HGVS/MolecularConsequence")
            if node.get("Type", "").strip()
        }
    )
    xrefs = archive.findall("./ClassifiedRecord/SimpleAllele/XRefList/XRef")
    rs_ids = sorted({f"rs{node.get('ID')}" for node in xrefs if node.get("DB") == "dbSNP" and node.get("Type") == "rs"})
    pubmed_ids = sorted(
        {
            node.text.strip()
            for node in classification.findall("./Citation/ID")
            if node.get("Source") == "PubMed" and node.text and node.text.strip()
        }
    )
    frequencies = []
    for node in archive.findall("./ClassifiedRecord/SimpleAllele/AlleleFrequencyList/AlleleFrequency"):
        try:
            frequencies.append(float(node.get("Value", "")))
        except ValueError:
            continue

    pathogenic = sum(class_counts.get(label, 0) for label in ("Pathogenic", "Likely pathogenic"))
    benign = sum(class_counts.get(label, 0) for label in ("Benign", "Likely benign"))
    vus = class_counts.get("Uncertain significance", 0)
    polarized = pathogenic > 0 and benign > 0
    submitters = int(classification.get("NumberOfSubmitters", archive.get("NumberOfSubmitters", "0")))
    review_status = (classification.findtext("ReviewStatus") or "").strip()
    criteria_provided = review_status.startswith("criteria provided")
    if criteria_provided and polarized and submitters >= 3:
        tier = "Tier 1: polarized, multi-submitter"
    elif criteria_provided and submitters >= 3:
        tier = "Tier 2: criteria-provided, multi-submitter"
    else:
        tier = "Tier 3: limited review support"

    accession = f"{archive.get('Accession')}.{archive.get('Version')}"
    variation_id = archive.get("VariationID", "")
    simple = archive.find("./ClassifiedRecord/SimpleAllele")
    return {
        "variation_id": variation_id,
        "accession_version": accession,
        "record_status": record_status,
        "variant_name": archive.get("VariationName", ""),
        "variant_type": archive.get("VariationType", ""),
        "canonical_spdi": (simple.findtext("CanonicalSPDI") if simple is not None else "") or "",
        "protein_change": (simple.findtext("ProteinChange") if simple is not None else "") or "",
        "rs_ids": "|".join(rs_ids),
        "molecular_consequences": "|".join(consequences),
        "grch38_accession": location.get("Accession", "") if location is not None else "",
        "grch38_position": location.get("positionVCF", "") if location is not None else "",
        "grch38_ref": location.get("referenceAlleleVCF", "") if location is not None else "",
        "grch38_alt": location.get("alternateAlleleVCF", "") if location is not None else "",
        "aggregate_classification": description,
        "review_status": review_status,
        "date_last_evaluated": date_text,
        "evaluation_age_years": round(age_years, 2) if not math.isnan(age_years) else "",
        "stale_over_five_years": bool(not math.isnan(age_years) and age_years > 5),
        "submission_count": int(classification.get("NumberOfSubmissions", archive.get("NumberOfSubmissions", "0"))),
        "submitter_count": submitters,
        "pathogenic_tier_count": pathogenic,
        "vus_count": vus,
        "benign_tier_count": benign,
        "other_classification_count": sum(class_counts.values()) - pathogenic - vus - benign,
        "classification_category_count": len(class_counts),
        "polarized_pathogenic_benign": polarized,
        "criteria_provided": criteria_provided,
        "priority_tier": tier,
        "trait_count": len(informative_traits),
        "traits": "|".join(informative_traits),
        "pubmed_citation_count": len(pubmed_ids),
        "pubmed_ids": "|".join(pubmed_ids),
        "maximum_reported_allele_frequency": max(frequencies) if frequencies else "",
        "clinvar_url": f"{CLINVAR_RECORD_BASE}/{variation_id}/",
        "snapshot_date": snapshot_date.isoformat(),
    }


def write_csv(path: Path, records: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not records:
        raise ValueError(f"Cannot write an empty table to {path}")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)


def build_classification_counts(records: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for record in records:
        for group, column in (
            ("Pathogenic or likely pathogenic", "pathogenic_tier_count"),
            ("Uncertain significance", "vus_count"),
            ("Benign or likely benign", "benign_tier_count"),
            ("Other", "other_classification_count"),
        ):
            rows.append(
                {
                    "accession_version": record["accession_version"],
                    "classification_group": group,
                    "submission_count": record[column],
                }
            )
    return rows


def summarize_records(records: list[dict[str, Any]]) -> dict[str, Any]:
    tiers = Counter(record["priority_tier"] for record in records)
    consequences = Counter(
        consequence
        for record in records
        for consequence in str(record["molecular_consequences"]).split("|")
        if consequence
    )
    return {
        "selected_records": len(records),
        "current_records": sum(record["record_status"] == "current" for record in records),
        "records_with_current_grch38_location": sum(bool(record["grch38_position"]) for record in records),
        "polarized_records": sum(record["polarized_pathogenic_benign"] for record in records),
        "stale_records": sum(record["stale_over_five_years"] for record in records),
        "criteria_provided_records": sum(record["criteria_provided"] for record in records),
        "median_submissions": _median([int(record["submission_count"]) for record in records]),
        "maximum_submissions": max(int(record["submission_count"]) for record in records),
        "tier_counts": dict(sorted(tiers.items())),
        "consequence_counts": dict(consequences.most_common()),
    }


def _median(values: list[int]) -> float:
    ordered = sorted(values)
    midpoint = len(ordered) // 2
    return float(ordered[midpoint]) if len(ordered) % 2 else (ordered[midpoint - 1] + ordered[midpoint]) / 2
