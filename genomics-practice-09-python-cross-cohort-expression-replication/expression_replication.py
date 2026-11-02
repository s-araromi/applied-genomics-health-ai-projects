"""Reusable functions for cross-platform paired-expression replication.

The functions operate on normalized values supplied by NCBI GEO. They preserve
sample-level pairing and never manufacture, impute, or silently substitute data.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import math
import re
import urllib.request
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


@dataclass(frozen=True)
class SourceFile:
    name: str
    url: str
    bytes: int
    md5: str


SOURCES = (
    SourceFile(
        "GSE19804_series_matrix.txt.gz",
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE19nnn/GSE19804/matrix/GSE19804_series_matrix.txt.gz",
        21_530_998,
        "c817eb7b8b74d1cbdd4db7580fb9ea08",
    ),
    SourceFile(
        "GPL570.annot.gz",
        "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL570/annot/GPL570.annot.gz",
        8_471_521,
        "be9fdbdb62e45257f2d3177a32f9a518",
    ),
    SourceFile(
        "GSE10072_series_matrix.txt.gz",
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE10nnn/GSE10072/matrix/GSE10072_series_matrix.txt.gz",
        13_021_909,
        "73273a4efa1d8e63777d7792b9bbabaf",
    ),
    SourceFile(
        "GPL96.annot.gz",
        "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL96/annot/GPL96.annot.gz",
        4_522_748,
        "a3c03c05e7b77d2305b84ba5e46181ec",
    ),
)


def md5sum(path: Path) -> str:
    digest = hashlib.md5(usedforsecurity=False)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def ensure_sources(cache: Path, *, allow_download: bool = True) -> None:
    """Download missing sources and reject any byte/hash mismatch."""
    cache.mkdir(parents=True, exist_ok=True)
    for source in SOURCES:
        path = cache / source.name
        if not path.exists():
            if not allow_download:
                raise FileNotFoundError(f"Required pinned source is absent: {path}")
            temporary = path.with_suffix(path.suffix + ".partial")
            request = urllib.request.Request(
                source.url,
                headers={"User-Agent": "applied-genomics-health-ai-projects/Project09"},
            )
            try:
                with urllib.request.urlopen(request, timeout=180) as response, temporary.open("wb") as output:
                    while block := response.read(1024 * 1024):
                        output.write(block)
                temporary.replace(path)
            except Exception:
                temporary.unlink(missing_ok=True)
                raise
        if path.stat().st_size != source.bytes or md5sum(path).lower() != source.md5:
            raise ValueError(f"Pinned source identity mismatch: {source.name}")


def _metadata_values(header: list[str], key: str) -> list[str]:
    lines = [line for line in header if line.startswith(key + "\t")]
    if len(lines) != 1:
        raise ValueError(f"Expected exactly one {key} metadata line")
    return next(csv.reader([lines[0].rstrip("\n")], delimiter="\t"))[1:]


def read_geo_matrix(
    path: Path,
    *,
    series: str,
    platform: str,
    expected_shape: tuple[int, int],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Parse one GEO Series Matrix and reconstruct participant/tissue labels."""
    header: list[str] = []
    with gzip.open(path, "rt", encoding="utf-8", errors="strict") as handle:
        for line in handle:
            if line.rstrip("\n") == "!series_matrix_table_begin":
                expression = pd.read_csv(handle, sep="\t", comment="!", index_col=0)
                break
            header.append(line)
        else:
            raise ValueError("GEO matrix table marker not found")
    if _metadata_values(header, "!Series_geo_accession") != [series]:
        raise ValueError("Series accession differs from expected source")
    if _metadata_values(header, "!Series_platform_id") != [platform]:
        raise ValueError("Platform accession differs from expected source")
    status = _metadata_values(header, "!Series_status")
    if len(status) != 1 or not status[0].startswith("Public on "):
        raise ValueError("GEO series is not marked public")
    titles = _metadata_values(header, "!Sample_title")
    accessions = _metadata_values(header, "!Sample_geo_accession")
    sources = _metadata_values(header, "!Sample_source_name_ch1")
    if not (len(titles) == len(accessions) == len(sources) == expression.shape[1]):
        raise ValueError("Sample metadata and expression columns differ in length")
    if expression.shape != expected_shape:
        raise ValueError(f"Unexpected {series} dimensions: {expression.shape}")
    if list(expression.columns) != accessions:
        raise ValueError("Expression columns are not aligned with GEO accessions")
    if expression.index.duplicated().any() or expression.columns.duplicated().any():
        raise ValueError("Duplicate probe or sample identifiers detected")
    values = expression.to_numpy(dtype=float)
    if not np.isfinite(values).all():
        raise ValueError("Missing or non-finite expression value detected")
    if values.min() < 0 or values.max() > 25:
        raise ValueError("Unexpected range for normalized log-scale microarray values")

    if series == "GSE19804":
        patient = pd.Series(titles).str.extract(r" ([0-9]+)[TN]$", expand=False)
        tissue = np.where(pd.Series(sources).str.contains("primary tumor", case=False), "Tumor",
                          np.where(pd.Series(sources).str.contains("adjacent normal", case=False), "Normal", ""))
    elif series == "GSE10072":
        patient = pd.Series(titles).str.extract(r"_(GT[0-9]+)$", expand=False)
        tissue = np.where(pd.Series(sources).eq("Adenocarcinoma of the Lung"), "Tumor",
                          np.where(pd.Series(sources).eq("Normal Lung Tissue"), "Normal", ""))
    else:
        raise ValueError(f"Unsupported series: {series}")
    metadata = pd.DataFrame(
        {"sample_accession": accessions, "sample_title": titles,
         "patient_id": patient, "tissue": tissue}
    )
    if metadata["patient_id"].isna().any() or (metadata["tissue"] == "").any():
        raise ValueError("Could not reconstruct participant and tissue labels")
    if metadata["sample_accession"].duplicated().any():
        raise ValueError("Duplicate sample accession detected")
    expression.index.name = "probe_id"
    return expression, metadata


def read_annotation(path: Path, expected_rows: int) -> pd.DataFrame:
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.rstrip("\n") == "!platform_table_begin":
                annotation = pd.read_csv(handle, sep="\t", comment="!", dtype=str, keep_default_na=False)
                break
        else:
            raise ValueError("GEO platform table marker not found")
    if len(annotation) != expected_rows or annotation["ID"].duplicated().any():
        raise ValueError("Unexpected platform annotation dimensions or duplicate probes")
    selected = annotation[["ID", "Gene ID", "Gene symbol", "Gene title"]].copy()
    selected.columns = ["probe_id", "gene_id", "gene_symbol", "gene_title"]
    eligible = selected["gene_id"].str.fullmatch(r"[0-9]+") & selected["gene_symbol"].str.fullmatch(r"[^/]+")
    return selected.loc[eligible].copy()


def representative_genes(expression: pd.DataFrame, annotation: pd.DataFrame) -> pd.DataFrame:
    """Select one probe per Entrez Gene ID by label-blind across-sample IQR."""
    mapped = annotation[annotation["probe_id"].isin(expression.index)].copy()
    if mapped.empty:
        raise ValueError("No annotated probes overlap the expression matrix")
    q75 = expression.loc[mapped["probe_id"]].quantile(0.75, axis=1)
    q25 = expression.loc[mapped["probe_id"]].quantile(0.25, axis=1)
    mapped["probe_iqr"] = mapped["probe_id"].map(q75 - q25)
    mapped = (mapped.sort_values(["gene_id", "probe_iqr", "probe_id"], ascending=[True, False, True])
              .drop_duplicates("gene_id", keep="first"))
    if mapped["gene_id"].duplicated().any():
        raise AssertionError("Representative gene table is not unique")
    return mapped.reset_index(drop=True)


def complete_pair_indices(metadata: pd.DataFrame) -> tuple[np.ndarray, np.ndarray, list[str]]:
    counts = metadata.groupby(["patient_id", "tissue"]).size().unstack(fill_value=0)
    eligible = sorted(counts.index[(counts.get("Tumor", 0) == 1) & (counts.get("Normal", 0) == 1)])
    tumor: list[int] = []
    normal: list[int] = []
    for patient in eligible:
        tumor.append(int(metadata.index[(metadata.patient_id == patient) & (metadata.tissue == "Tumor")][0]))
        normal.append(int(metadata.index[(metadata.patient_id == patient) & (metadata.tissue == "Normal")][0]))
    return np.asarray(tumor), np.asarray(normal), eligible


def benjamini_hochberg(p_values: np.ndarray) -> np.ndarray:
    p = np.asarray(p_values, dtype=float)
    if p.ndim != 1 or np.any(~np.isfinite(p)) or np.any((p < 0) | (p > 1)):
        raise ValueError("p-values must be a finite one-dimensional array in [0, 1]")
    order = np.argsort(p, kind="stable")
    ranked = p[order] * len(p) / np.arange(1, len(p) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adjusted = np.empty_like(ranked)
    adjusted[order] = np.clip(ranked, 0, 1)
    return adjusted


def paired_statistics(values: np.ndarray, tumor: np.ndarray, normal: np.ndarray) -> pd.DataFrame:
    """Vectorized paired mean difference, t test, standard error and CI."""
    if values.ndim != 2 or len(tumor) != len(normal) or len(tumor) < 2:
        raise ValueError("Expected genes by samples with at least two aligned pairs")
    difference = values[:, tumor] - values[:, normal]
    n = difference.shape[1]
    effect = difference.mean(axis=1)
    sd = difference.std(axis=1, ddof=1)
    se = sd / math.sqrt(n)
    statistic = np.divide(effect, se, out=np.zeros_like(effect), where=se > 0)
    p_value = 2 * stats.t.sf(np.abs(statistic), df=n - 1)
    zero = se == 0
    p_value[zero & (effect != 0)] = 0
    p_value[zero & (effect == 0)] = 1
    critical = stats.t.ppf(0.975, df=n - 1)
    return pd.DataFrame({
        "n_pairs": n, "mean_log2_difference": effect, "standard_error": se,
        "ci95_lower": effect - critical * se, "ci95_upper": effect + critical * se,
        "t_statistic": statistic, "p_value": p_value,
    })


def clustered_ols(values: np.ndarray, metadata: pd.DataFrame) -> pd.DataFrame:
    """All-sample tissue contrast with participant-clustered sandwich SEs.

    The small-sample correction is G/(G-1) * (N-1)/(N-k), with t reference
    degrees of freedom G-1. Each participant contributes one cluster.
    """
    y = np.asarray(values, dtype=float).T  # samples by genes
    if y.shape[0] != len(metadata) or not np.isfinite(y).all():
        raise ValueError("Values and sample metadata are not aligned")
    tissue = (metadata["tissue"].to_numpy() == "Tumor").astype(float)
    x = np.column_stack((np.ones(len(tissue)), tissue))
    inv = np.linalg.inv(x.T @ x)
    beta = inv @ x.T @ y
    residual = y - x @ beta
    v = inv[:, 1]
    score_squares = np.zeros(y.shape[1])
    clusters = metadata["patient_id"].astype(str).to_numpy()
    unique = np.unique(clusters)
    for cluster in unique:
        take = clusters == cluster
        # v' X_g' u_g = (X_g v)' u_g, evaluated for every gene.
        score = (x[take] @ v) @ residual[take]
        score_squares += score**2
    n, k, groups = len(metadata), x.shape[1], len(unique)
    correction = groups / (groups - 1) * (n - 1) / (n - k)
    se = np.sqrt(np.maximum(0, correction * score_squares))
    statistic = np.divide(beta[1], se, out=np.zeros_like(se), where=se > 0)
    p_value = 2 * stats.t.sf(np.abs(statistic), df=groups - 1)
    zero = se == 0
    p_value[zero & (beta[1] != 0)] = 0
    p_value[zero & (beta[1] == 0)] = 1
    return pd.DataFrame({"clustered_effect": beta[1], "clustered_standard_error": se,
                         "clustered_t": statistic, "clustered_p_value": p_value,
                         "n_samples": n, "n_participants": groups})


def fixed_effect_meta(effect_1: np.ndarray, se_1: np.ndarray,
                      effect_2: np.ndarray, se_2: np.ndarray) -> pd.DataFrame:
    arrays = [np.asarray(x, dtype=float) for x in (effect_1, se_1, effect_2, se_2)]
    if any(x.shape != arrays[0].shape for x in arrays) or np.any(arrays[1] <= 0) or np.any(arrays[3] <= 0):
        raise ValueError("Meta-analysis inputs must align and have positive standard errors")
    w1, w2 = 1 / arrays[1] ** 2, 1 / arrays[3] ** 2
    effect = (w1 * arrays[0] + w2 * arrays[2]) / (w1 + w2)
    se = np.sqrt(1 / (w1 + w2))
    z = effect / se
    p_value = 2 * stats.norm.sf(np.abs(z))
    q = w1 * (arrays[0] - effect) ** 2 + w2 * (arrays[2] - effect) ** 2
    i2 = np.where(q > 0, np.maximum(0, (q - 1) / q), 0)
    return pd.DataFrame({"fixed_effect_log2_difference": effect,
                         "fixed_effect_standard_error": se,
                         "fixed_effect_z": z, "fixed_effect_p_value": p_value,
                         "heterogeneity_Q": q, "heterogeneity_I2": i2})


def validate_output_accounting(results: pd.DataFrame) -> None:
    if results["gene_id"].duplicated().any() or results.isna().any().any():
        raise ValueError("Result table contains duplicate genes or missing values")
    if not np.all((results["discovery_q"] >= 0) & (results["discovery_q"] <= 1)):
        raise ValueError("Invalid discovery q-values")
    if (results["replicated"] & ~results["discovery_candidate"]).any():
        raise ValueError("Replicated genes must be discovery candidates")
