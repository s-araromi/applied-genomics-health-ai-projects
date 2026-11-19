#!/usr/bin/env python3
"""Retrieve pinned GSE2034 records and regenerate the checked analytical subset."""
from __future__ import annotations

import argparse
import shutil
import tempfile
import urllib.request
from pathlib import Path

import pandas as pd

from breast_relapse_ml import (
    collapse_to_unambiguous_genes,
    construct_five_year_endpoint,
    parse_clinical_table,
    parse_platform_annotation,
    parse_series_matrix,
    sha256,
)

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
SOURCES = {
    "GSE2034_series_matrix.txt.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE2nnn/GSE2034/matrix/GSE2034_series_matrix.txt.gz",
        "f6a99c62f0a291858f14806ad7b5faf1da159f88630d11cfe994815ff3bf9ea3",
    ),
    "GSE2034_family.soft.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE2nnn/GSE2034/soft/GSE2034_family.soft.gz",
        "1fa927dc72870434ba799041dcb6dcd292cbb77ab1eac203ea8e5d51272cd279",
    ),
    "GPL96.annot.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL96/annot/GPL96.annot.gz",
        "88e0b22362bac779eb220b3b185c80faa6510a92b9358eaad159a561ab4351c4",
    ),
}


def download(url: str, target: Path) -> None:
    request = urllib.request.Request(url, headers={"User-Agent": "Project-11-reproducibility/1.0"})
    with urllib.request.urlopen(request, timeout=120) as response, target.open("wb") as output:
        shutil.copyfileobj(response, output)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, default=DATA / "cache")
    parser.add_argument("--no-download", action="store_true")
    args = parser.parse_args()
    args.cache.mkdir(parents=True, exist_ok=True)
    for name, (url, expected) in SOURCES.items():
        path = args.cache / name
        if not path.exists():
            if args.no_download:
                raise FileNotFoundError(f"Pinned cache file is absent: {path}")
            with tempfile.NamedTemporaryFile(dir=args.cache, delete=False) as tmp:
                temporary = Path(tmp.name)
            try:
                download(url, temporary)
                temporary.replace(path)
            finally:
                temporary.unlink(missing_ok=True)
        observed = sha256(path)
        if observed != expected:
            raise ValueError(f"Source identity mismatch for {name}: {observed}")

    matrix, metadata = parse_series_matrix(args.cache / "GSE2034_series_matrix.txt.gz")
    clinical = construct_five_year_endpoint(parse_clinical_table(args.cache / "GSE2034_family.soft.gz"))
    annotation = parse_platform_annotation(args.cache / "GPL96.annot.gz")
    genes, selected = collapse_to_unambiguous_genes(matrix, annotation)
    accessions = metadata.get("!Sample_geo_accession", [])
    if list(matrix.columns) != accessions or set(clinical["gsm"]) != set(matrix.columns):
        raise ValueError("Expression and clinical sample identities do not reconcile")

    clinical_out = clinical.rename(columns={"PID": "participant_id"})
    clinical_out.to_csv(DATA / "clinical_metadata.csv", index=False, lineterminator="\n")
    selected.to_csv(DATA / "representative_probes.csv", index=False, lineterminator="\n")
    genes.to_csv(DATA / "expression_genes.csv.gz", compression={"method": "gzip", "mtime": 0},
                 lineterminator="\n")
    print(f"Verified 286 samples, {matrix.shape[0]:,} probes, and {genes.shape[0]:,} genes")
    print("Regenerated clinical_metadata.csv, representative_probes.csv, and expression_genes.csv.gz")


if __name__ == "__main__":
    main()
