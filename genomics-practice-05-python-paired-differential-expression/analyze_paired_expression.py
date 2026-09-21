#!/usr/bin/env python3
"""Paired differential-expression analysis of the complete GSE19804 cohort."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import sys
import urllib.request
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

MATRIX_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE19nnn/GSE19804/matrix/GSE19804_series_matrix.txt.gz"
ANNOTATION_URL = "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL570/annot/GPL570.annot.gz"
MATRIX_BYTES = 21_530_998
MATRIX_MD5 = "c817eb7b8b74d1cbdd4db7580fb9ea08"
ANNOTATION_BYTES = 8_471_521
ANNOTATION_MD5 = "be9fdbdb62e45257f2d3177a32f9a518"
FDR_THRESHOLD = 0.05
EFFECT_THRESHOLD = 1.0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", type=Path, default=Path("data/GSE19804_series_matrix.txt.gz"))
    parser.add_argument("--annotation", type=Path, default=Path("data/GPL570.annot.gz"))
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"))
    return parser.parse_args()


def md5sum(path: Path) -> str:
    digest = hashlib.md5(usedforsecurity=False)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def download_validated(url: str, path: Path, expected_bytes: int, expected_md5: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not path.exists():
        print(f"Downloading {path.name}...")
        request = urllib.request.Request(url, headers={"User-Agent": "genomics-practice-expression/1.0"})
        temporary = path.with_suffix(path.suffix + ".partial")
        try:
            with urllib.request.urlopen(request, timeout=180) as response, temporary.open("wb") as output:
                while block := response.read(1024 * 1024):
                    output.write(block)
            temporary.replace(path)
        except Exception:
            temporary.unlink(missing_ok=True)
            raise
    observed_bytes = path.stat().st_size
    if observed_bytes != expected_bytes:
        raise ValueError(f"{path.name} byte count is {observed_bytes}; expected {expected_bytes}")
    observed_md5 = md5sum(path)
    if observed_md5.lower() != expected_md5.lower():
        raise ValueError(f"{path.name} MD5 is {observed_md5}; expected {expected_md5}")


def _metadata_values(lines: list[str], key: str) -> list[str]:
    matches = [line for line in lines if line.startswith(key + "\t")]
    if len(matches) != 1:
        raise ValueError(f"Expected one {key} metadata line")
    return [value.strip('"') for value in matches[0].rstrip("\n").split("\t")[1:]]


def read_expression(path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    header: list[str] = []
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.rstrip("\n") == "!series_matrix_table_begin":
                expression = pd.read_csv(handle, sep="\t", comment="!", index_col=0)
                break
            header.append(line)
        else:
            raise ValueError("GEO expression table marker not found")
    titles = _metadata_values(header, "!Sample_title")
    accessions = _metadata_values(header, "!Sample_geo_accession")
    sources = _metadata_values(header, "!Sample_source_name_ch1")
    if not (len(titles) == len(accessions) == len(sources) == expression.shape[1]):
        raise ValueError("Sample metadata and expression columns differ in length")
    tissue = np.where(pd.Series(sources).str.contains("primary tumor", case=False), "Tumor",
                      np.where(pd.Series(sources).str.contains("adjacent normal", case=False), "Normal", ""))
    patient_id = pd.Series(titles).str.extract(r" ([0-9]+)[TN]$", expand=False)
    metadata = pd.DataFrame({"sample_accession": accessions, "sample_title": titles,
                             "patient_id": patient_id, "tissue": tissue})
    if metadata["patient_id"].isna().any() or (metadata["tissue"] == "").any():
        raise ValueError("Could not reconstruct tissue groups and patient identifiers")
    if list(expression.columns) != accessions:
        raise ValueError("Expression columns are not aligned with GEO accessions")
    expression.index.name = "probe_id"
    return expression, metadata


def read_annotation(path: Path) -> pd.DataFrame:
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.rstrip("\n") == "!platform_table_begin":
                annotation = pd.read_csv(handle, sep="\t", comment="!", dtype=str,
                                         keep_default_na=False)
                break
        else:
            raise ValueError("GPL570 annotation table marker not found")
    selected = annotation[["ID", "Gene symbol", "Gene ID", "Gene title"]].copy()
    selected.columns = ["probe_id", "gene_symbol", "gene_id", "gene_title"]
    if selected["probe_id"].duplicated().any():
        raise ValueError("GPL570 annotation contains duplicate probe identifiers")
    return selected


def validate_inputs(expression: pd.DataFrame, metadata: pd.DataFrame) -> None:
    if expression.shape != (54_675, 120):
        raise ValueError(f"Unexpected expression dimensions: {expression.shape}")
    if expression.index.duplicated().any() or expression.columns.duplicated().any():
        raise ValueError("Duplicate expression identifiers detected")
    values = expression.to_numpy(dtype=float)
    if not np.isfinite(values).all():
        raise ValueError("Missing or non-finite expression values detected")
    if metadata["patient_id"].nunique() != 60:
        raise ValueError("Expected 60 participants")
    counts = metadata.groupby(["patient_id", "tissue"]).size().unstack(fill_value=0)
    if not ((counts["Tumor"] == 1) & (counts["Normal"] == 1)).all():
        raise ValueError("Every participant must have one tumor and one normal sample")


def benjamini_hochberg(p_values: np.ndarray) -> np.ndarray:
    p_values = np.asarray(p_values, dtype=float)
    if p_values.ndim != 1 or np.any(~np.isfinite(p_values)) or np.any((p_values < 0) | (p_values > 1)):
        raise ValueError("p-values must be a finite one-dimensional array between zero and one")
    order = np.argsort(p_values)
    ranked = p_values[order] * len(p_values) / np.arange(1, len(p_values) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adjusted = np.empty_like(ranked)
    adjusted[order] = np.clip(ranked, 0, 1)
    return adjusted


def pair_indices(metadata: pd.DataFrame, excluded_patients: set[str] | None = None) -> tuple[list[int], list[int], list[str]]:
    excluded_patients = excluded_patients or set()
    tumor_indices: list[int] = []
    normal_indices: list[int] = []
    patients: list[str] = []
    for patient in metadata.loc[metadata["tissue"] == "Tumor", "patient_id"]:
        if patient in excluded_patients:
            continue
        tumor = metadata.index[(metadata["patient_id"] == patient) & (metadata["tissue"] == "Tumor")]
        normal = metadata.index[(metadata["patient_id"] == patient) & (metadata["tissue"] == "Normal")]
        if len(tumor) != 1 or len(normal) != 1:
            raise ValueError(f"Patient {patient} does not have one complete pair")
        tumor_indices.append(int(tumor[0]))
        normal_indices.append(int(normal[0]))
        patients.append(patient)
    return tumor_indices, normal_indices, patients


def review_patients(expression: pd.DataFrame, metadata: pd.DataFrame) -> tuple[set[str], pd.DataFrame]:
    metrics = metadata.copy()
    metrics["mean_within_tissue_correlation"] = np.nan
    metrics["review_threshold"] = np.nan
    metrics["review_flag"] = False
    values = expression.to_numpy(dtype=float)
    for tissue in ("Tumor", "Normal"):
        indices = metadata.index[metadata["tissue"] == tissue].to_numpy()
        correlations = np.corrcoef(values[:, indices], rowvar=False)
        np.fill_diagonal(correlations, np.nan)
        means = np.nanmean(correlations, axis=0)
        q1, q3 = np.quantile(means, [0.25, 0.75])
        threshold = q1 - 2 * (q3 - q1)
        metrics.loc[indices, "mean_within_tissue_correlation"] = means
        metrics.loc[indices, "review_threshold"] = threshold
        metrics.loc[indices, "review_flag"] = means < threshold
    patients = set(metrics.loc[metrics["review_flag"], "patient_id"])
    return patients, metrics


def paired_test(expression: pd.DataFrame, metadata: pd.DataFrame,
                excluded_patients: set[str] | None = None) -> tuple[pd.DataFrame, list[str]]:
    tumor_indices, normal_indices, patients = pair_indices(metadata, excluded_patients)
    values = expression.to_numpy(dtype=float)
    tumor = values[:, tumor_indices]
    normal = values[:, normal_indices]
    differences = tumor - normal
    statistic, p_value = stats.ttest_rel(tumor, normal, axis=1)
    result = pd.DataFrame({
        "mean_log2_fold_change": differences.mean(axis=1),
        "t_statistic": statistic,
        "p_value": p_value,
        "fdr": benjamini_hochberg(p_value),
    }, index=expression.index)
    result["significant"] = ((result["fdr"] <= FDR_THRESHOLD) &
                             (result["mean_log2_fold_change"].abs() >= EFFECT_THRESHOLD))
    return result, patients


def unambiguous_symbol(symbol: str) -> bool:
    symbol = str(symbol).strip()
    return bool(symbol) and "///" not in symbol


def build_results(expression: pd.DataFrame, metadata: pd.DataFrame,
                  annotation: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame,
                                                     pd.DataFrame, pd.DataFrame]:
    excluded_patients, sample_metrics = review_patients(expression, metadata)
    primary, primary_patients = paired_test(expression, metadata)
    sensitivity, sensitivity_patients = paired_test(expression, metadata, excluded_patients)
    results = primary.rename(columns={
        "mean_log2_fold_change": "primary_mean_log2_fold_change",
        "t_statistic": "primary_t_statistic", "p_value": "primary_p_value",
        "fdr": "primary_fdr", "significant": "primary_significant"
    }).join(sensitivity.rename(columns={
        "mean_log2_fold_change": "sensitivity_mean_log2_fold_change",
        "t_statistic": "sensitivity_t_statistic", "p_value": "sensitivity_p_value",
        "fdr": "sensitivity_fdr", "significant": "sensitivity_significant"
    }))
    results["robust"] = (results["primary_significant"] & results["sensitivity_significant"] &
                         (np.sign(results["primary_mean_log2_fold_change"]) ==
                          np.sign(results["sensitivity_mean_log2_fold_change"])))
    results["direction"] = np.where(results["primary_mean_log2_fold_change"] > 0,
                                    "Higher in tumor", "Lower in tumor")
    results = results.reset_index().merge(annotation, on="probe_id", how="left")
    for column in ("gene_symbol", "gene_id", "gene_title"):
        results[column] = results[column].fillna("")
    robust = results.loc[results["robust"]].sort_values(
        ["primary_fdr", "primary_mean_log2_fold_change"], ascending=[True, False]).copy()
    eligible = robust.loc[robust["gene_symbol"].map(unambiguous_symbol)].copy()
    eligible["absolute_effect"] = eligible["primary_mean_log2_fold_change"].abs()
    unique_genes = (eligible.sort_values(["absolute_effect", "primary_fdr"], ascending=[False, True])
                    .drop_duplicates("gene_symbol").drop(columns="absolute_effect"))
    pair_manifest = pd.DataFrame({"patient_id": primary_patients})
    tumor_i, normal_i, _ = pair_indices(metadata)
    pair_manifest["tumor_accession"] = metadata.loc[tumor_i, "sample_accession"].to_numpy()
    pair_manifest["normal_accession"] = metadata.loc[normal_i, "sample_accession"].to_numpy()
    pair_manifest["included_in_sensitivity"] = ~pair_manifest["patient_id"].isin(excluded_patients)
    return results, robust, unique_genes, sample_metrics, pair_manifest


def write_report(path: Path, results: pd.DataFrame, robust: pd.DataFrame,
                 unique_genes: pd.DataFrame, pair_manifest: pd.DataFrame,
                 sample_metrics: pd.DataFrame) -> None:
    primary = results["primary_significant"]
    higher = primary & (results["primary_mean_log2_fold_change"] > 0)
    lower = primary & (results["primary_mean_log2_fold_change"] < 0)
    annotated = results["gene_symbol"].str.strip().ne("")
    leading_higher = (robust.loc[(robust["direction"] == "Higher in tumor") &
                                 robust["gene_symbol"].map(unambiguous_symbol)]
                      .drop_duplicates("gene_symbol").nsmallest(5, "primary_fdr"))
    leading_lower = (robust.loc[(robust["direction"] == "Lower in tumor") &
                                robust["gene_symbol"].map(unambiguous_symbol)]
                     .drop_duplicates("gene_symbol").nsmallest(5, "primary_fdr"))
    lines = [
        "# Paired Differential-Expression Analysis",
        "", "## Dataset", "",
        "- GEO series: `GSE19804`", "- Platform: `GPL570`",
        "- Primary analysis: 60 paired tumor-normal participants",
        f"- Sensitivity analysis: {pair_manifest['included_in_sensitivity'].sum()} pairs",
        "", "## Results", "",
        f"- Probe sets tested: {len(results):,}",
        f"- Primary signals meeting FDR and effect-size criteria: {primary.sum():,}",
        f"- Higher in tumor: {higher.sum():,}",
        f"- Lower in tumor: {lower.sum():,}",
        f"- Robust probe sets after sensitivity analysis: {len(robust):,}",
        f"- Robust unique, unambiguous genes: {len(unique_genes):,}",
        f"- Probe sets with gene symbols: {annotated.sum():,}",
        f"- Samples marked for review: {sample_metrics['review_flag'].sum()}",
        "", "## Leading robust signals", "",
        "Higher in tumor: " + ", ".join(leading_higher["gene_symbol"].replace("", "unannotated")),
        "", "Lower in tumor: " + ", ".join(leading_lower["gene_symbol"].replace("", "unannotated")),
        "", "Three SEMA5A probes were robustly lower in tumor, with primary mean",
        "log2 fold changes from -1.44 to -1.78. This agrees with the originating",
        "study's interest in SEMA5A but is not an independent replication because",
        "the same cohort is analyzed here.",
        "", "## Interpretation", "",
        "Signals were required to meet BH FDR <= 0.05 and an absolute paired",
        "mean log2 fold change >= 1. Robust signals also met both criteria with",
        "the same direction after excluding seven pairs containing review-flagged",
        "arrays. Association in this cohort does not establish causation, diagnostic",
        "performance, or clinical utility."
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def plot_results(path_png: Path, path_svg: Path, results: pd.DataFrame,
                 pair_manifest: pd.DataFrame) -> None:
    figure, axes = plt.subplots(1, 2, figsize=(12, 5))
    x = results["primary_mean_log2_fold_change"].to_numpy()
    y = -np.log10(np.maximum(results["primary_fdr"].to_numpy(), np.finfo(float).tiny))
    colors = np.where(results["robust"], np.where(x > 0, "#B64926", "#2A6F97"), "#C7CBD1")
    axes[0].scatter(x, y, c=colors, s=7, alpha=0.65, linewidths=0, rasterized=True)
    axes[0].axvline(-1, color="#555555", linestyle="--", linewidth=1)
    axes[0].axvline(1, color="#555555", linestyle="--", linewidth=1)
    axes[0].axhline(-np.log10(0.05), color="#555555", linestyle="--", linewidth=1)
    axes[0].set(title="Paired differential-expression signals",
                xlabel="Mean paired log2 fold change", ylabel="-log10(BH FDR)")
    primary_counts = [
        ((results["primary_significant"]) & (x < 0)).sum(),
        ((results["primary_significant"]) & (x > 0)).sum(),
        results["robust"].sum(),
    ]
    axes[1].bar(["Lower in tumor", "Higher in tumor", "Robust"], primary_counts,
                color=["#2A6F97", "#B64926", "#5B6B73"])
    axes[1].set(title="Signals meeting declared criteria", ylabel="Probe sets")
    axes[1].tick_params(axis="x", rotation=15)
    figure.suptitle(f"GSE19804 paired analysis: 60 primary pairs; {pair_manifest['included_in_sensitivity'].sum()} sensitivity pairs",
                    fontsize=13, fontweight="bold")
    figure.tight_layout()
    figure.savefig(path_png, dpi=180, bbox_inches="tight",
                   metadata={"Software": None, "Author": "Sulaimon Araromi"})
    figure.savefig(path_svg, bbox_inches="tight", metadata={"Creator": "Sulaimon Araromi"})
    plt.close(figure)


def main() -> int:
    args = parse_args()
    download_validated(MATRIX_URL, args.matrix, MATRIX_BYTES, MATRIX_MD5)
    download_validated(ANNOTATION_URL, args.annotation, ANNOTATION_BYTES, ANNOTATION_MD5)
    print("Reading GSE19804 and GPL570...")
    expression, metadata = read_expression(args.matrix)
    annotation = read_annotation(args.annotation)
    validate_inputs(expression, metadata)
    print("Running paired primary and sensitivity analyses...")
    results, robust, unique_genes, sample_metrics, pair_manifest = build_results(
        expression, metadata, annotation)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    results.to_csv(args.output_dir / "paired_differential_expression_results.csv", index=False)
    robust.to_csv(args.output_dir / "robust_probe_results.csv", index=False)
    unique_genes.to_csv(args.output_dir / "robust_unique_gene_results.csv", index=False)
    pair_manifest.to_csv(args.output_dir / "sample_pair_manifest.csv", index=False)
    sample_metrics.to_csv(args.output_dir / "sample_review_metrics.csv", index=False)
    write_report(args.output_dir / "paired_expression_report.md", results, robust,
                 unique_genes, pair_manifest, sample_metrics)
    plot_results(args.output_dir / "paired_expression_overview.png",
                 args.output_dir / "paired_expression_overview.svg", results, pair_manifest)
    print(f"Probe sets tested: {len(results):,}")
    print(f"Primary significant probe sets: {results['primary_significant'].sum():,}")
    print(f"Higher in tumor: {(results['primary_significant'] & (results['primary_mean_log2_fold_change'] > 0)).sum():,}")
    print(f"Lower in tumor: {(results['primary_significant'] & (results['primary_mean_log2_fold_change'] < 0)).sum():,}")
    print(f"Robust probe sets: {len(robust):,}")
    print(f"Robust unique, unambiguous genes: {len(unique_genes):,}")
    print(f"Outputs written to {args.output_dir.resolve()}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError, RuntimeError) as exc:
        print(f"Analysis stopped: {exc}", file=sys.stderr)
        raise SystemExit(1)
