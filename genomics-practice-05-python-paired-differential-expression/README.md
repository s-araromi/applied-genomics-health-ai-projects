# Project 05: Paired Differential-Expression Analysis

## Project summary

This Python project extends Project 04 by testing expression differences between
60 paired lung-tumor and adjacent-normal samples from
[`GSE19804`](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE19804). It
integrates the complete processed matrix with the official
[`GPL570`](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL570) annotation.

The workflow uses participant-matched tests, Benjamini–Hochberg false-discovery
rate control, an explicit effect-size threshold, and a sensitivity analysis that
omits pairs containing arrays marked for review in Project 04.

## Analytical objectives

- Verify the exact GEO matrix and GPL570 annotation with byte counts and MD5.
- Reconstruct all 60 participant-matched tumor-normal pairs.
- Perform a paired t-test for every probe set.
- Control the false-discovery rate across 54,675 tests.
- Require both statistical evidence and a biologically interpretable effect size.
- Re-run the analysis after excluding review-flagged pairs.
- Integrate gene symbols, identifiers, and titles from GPL570.
- Produce probe-level, robust, and unique-gene result tables.

## Data provenance

| File | Bytes | MD5 |
| --- | ---: | --- |
| `GSE19804_series_matrix.txt.gz` | 21,530,998 | `c817eb7b8b74d1cbdd4db7580fb9ea08` |
| `GPL570.annot.gz` | 8,471,521 | `be9fdbdb62e45257f2d3177a32f9a518` |

Both official files are downloaded from NCBI GEO at runtime and excluded from
Git. Reuse is subject to the
[NCBI data-usage policies](https://www.ncbi.nlm.nih.gov/home/about/policies/).

## Statistical criteria

A primary signal must meet both:

- Benjamini–Hochberg FDR ≤ 0.05.
- Absolute paired mean log2 fold change ≥ 1.

A robust signal must meet the same criteria and direction in the 53-pair
sensitivity analysis. Seven pairs are omitted from that analysis because at
least one array was marked for review in Project 04; all 60 pairs remain in the
primary analysis.

## Verified results

| Result | Value |
| --- | ---: |
| Probe sets tested | 54,675 |
| Primary tumor-normal pairs | 60 |
| Sensitivity-analysis pairs | 53 |
| Primary significant probe sets | 1,960 |
| Higher in tumor | 591 |
| Lower in tumor | 1,369 |
| Robust probe sets | 1,928 |
| Robust unique, unambiguous genes | 1,328 |
| Probe sets with gene symbols | 45,118 |

Leading robust signals included lower tumor expression for probes annotated to
`SPTBN1`, `AGER`, `CA4`, `SPOCK2`, and `GPM6A`, and higher expression for probes
annotated to `GOLM1`, `OCIAD2`, `PYCR1`, `COL10A1`, and `SPP1`.

Three `SEMA5A` probes were robustly lower in tumor, with mean log2 fold changes
between −1.44 and −1.78. This agrees with the originating study's interest in
`SEMA5A`, but it is not an independent replication because the same cohort is
analyzed here.

## Requirements

- Python 3.10 or newer.
- `numpy`, `pandas`, `scipy`, and `matplotlib`.
- Internet access for the first source-file download.

## Reproduction

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
python analyze_paired_expression.py
python -m unittest discover -s tests -v
```

## Outputs

| File | Description |
| --- | --- |
| `outputs/paired_differential_expression_results.csv` | Complete 54,675-row statistical result |
| `outputs/robust_probe_results.csv` | Probe sets robust to sensitivity analysis |
| `outputs/robust_unique_gene_results.csv` | One representative per unambiguous gene symbol |
| `outputs/sample_pair_manifest.csv` | Pair membership and sensitivity inclusion |
| `outputs/sample_review_metrics.csv` | Reconstructed Project 04 correlation flags |
| `outputs/paired_expression_report.md` | Findings and interpretation |
| `outputs/paired_expression_overview.png` | Raster statistical visualization |
| `outputs/paired_expression_overview.svg` | Vector statistical visualization |

## Interpretation limits

Differential expression is an association within this cohort. It does not prove
causation, diagnostic accuracy, therapeutic response, or clinical utility.
Multiple probe sets can represent one gene, and platform annotation can contain
ambiguous mappings.
