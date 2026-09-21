# Project 04: Paired GEO Expression-Matrix Quality Assessment

## Project summary

This R project evaluates the complete processed expression matrix from
[`GSE19804`](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE19804), a
paired lung-cancer study of 60 never-smoking women. Each participant contributed
a primary lung-tumor sample and adjacent-normal lung tissue, producing 120 arrays
on the [`GPL570`](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL570)
platform.

The project focuses on data integrity, distributional quality, sample similarity,
and paired-study structure before differential-expression analysis.

## Analytical objectives

- Download and checksum-validate the complete GEO series matrix.
- Confirm matrix dimensions, identifiers, numerical completeness, and value range.
- Reconstruct all 60 tumor-normal participant pairs from official metadata.
- Compare sample medians and within-sample expression spread.
- Calculate within-tissue and within-pair Pearson correlations.
- Mark unusual samples for review without automatically excluding them.
- Produce reproducible tables, a report, and publication-ready graphics.

## Data provenance

| Field | Value |
| --- | --- |
| GEO series | `GSE19804` |
| Platform | `GPL570` |
| Matrix | `GSE19804_series_matrix.txt.gz` |
| Compressed size | 21,530,998 bytes |
| MD5 | `c817eb7b8b74d1cbdd4db7580fb9ea08` |
| Participants | 60 |
| Samples | 120 |
| Probe sets | 54,675 |

The complete matrix is downloaded when the analysis runs and is excluded from
Git. Data use remains subject to the
[NCBI data-usage policies](https://www.ncbi.nlm.nih.gov/home/about/policies/).

## Verified results

| Assessment | Result |
| --- | ---: |
| Probe sets | 54,675 |
| Samples | 120 |
| Complete tumor-normal pairs | 60 |
| Missing or non-finite values | 0 |
| Duplicate probe identifiers | 0 |
| Duplicate sample identifiers | 0 |
| Expression range | 3.0427–14.8903 |
| Sample-median range | 6.02783–6.02919 |
| Median within-pair correlation | 0.9489 |
| Samples marked for review | 8 |
| Automatically excluded | 0 |

The narrow sample-median range is consistent with the quantile normalization
reported in the GEO matrix. Eight arrays had unusually low mean correlation with
other arrays from the same tissue group under the declared Q1 minus 2 × IQR
rule. They remain in the dataset because technical variation, tumor
heterogeneity, and participant-specific biology cannot be distinguished by this
QC screen alone.

## Requirements

- R 4.3 or newer.
- Internet access for the first matrix download.
- No external R packages.

## Reproduction

Run the full analysis:

```bash
Rscript analyze_geo_expression_qc.R
```

The validated matrix is saved in `data/` and reused on subsequent runs.

Run the offline automated checks:

```bash
Rscript tests/test_geo_expression_qc.R
```

An existing validated matrix can be supplied explicitly:

```bash
Rscript analyze_geo_expression_qc.R --input path/to/GSE19804_series_matrix.txt.gz
```

## Outputs

| File | Description |
| --- | --- |
| `outputs/expression_qc_summary.csv` | Dataset-level integrity and QC statistics |
| `outputs/sample_quality_metrics.csv` | Per-array distributions, correlations, and flags |
| `outputs/pair_correlations.csv` | Tumor-normal correlations for 60 participants |
| `outputs/expression_qc_report.md` | Results and interpretation |
| `outputs/expression_qc_overview.png` | Raster QC visualization |
| `outputs/expression_qc_overview.svg` | Vector QC visualization |

## Interpretation limits

Review flags are exploratory and are not automatic exclusion decisions. This
project does not test differential expression, diagnose lung cancer, infer
causality, or establish clinical biomarkers.
