# Analysis Workflow: Paired GEO Expression-Matrix Quality Assessment

## 1. Objective

Evaluate the integrity and sample-level quality of the complete processed
`GSE19804` expression matrix while preserving its paired tumor-normal design.

## 2. Source acquisition and verification

The workflow downloads the official compressed series matrix from NCBI GEO. It
accepts the file only when both the pinned byte count and MD5 checksum match:

```text
Bytes: 21,530,998
MD5: c817eb7b8b74d1cbdd4db7580fb9ea08
```

The download is cached locally under `data/` and excluded from version control.

## 3. Metadata reconstruction

The program extracts sample titles, GEO accessions, and tissue-source fields from
the matrix header. Titles such as `Lung Cancer 103T` and `Lung Normal 103N` are
used to recover participant identifiers and verify exactly one tumor and one
normal sample per participant.

## 4. Matrix validation

The processed expression table is required to contain:

- 54,675 unique probe-set rows.
- 120 unique GEO sample columns.
- Sample columns aligned with matrix-header accessions.
- Numeric and finite expression values throughout.
- 60 tumor samples and 60 adjacent-normal samples.

## 5. Distributional assessment

For every array, the workflow calculates the median, interquartile range,
minimum, and maximum expression value. The observed sample medians are expected
to be very similar because GEO reports that the processed data were quantile
normalized.

## 6. Correlation assessment

Two complementary correlations are calculated:

- Mean within-tissue correlation: each sample's mean Pearson correlation with
  the other 59 arrays from the same tissue group.
- Within-pair correlation: Pearson correlation between tumor and adjacent-normal
  expression profiles from the same participant.

## 7. Review flags

Within each tissue group, the lower review boundary is:

```text
Q1 − 2 × IQR
```

Samples below this boundary are marked for review. No samples are automatically
removed. This preserves biological heterogeneity and keeps downstream sensitivity
analysis possible.

## 8. Outputs

The workflow writes three CSV files, one Markdown report, and equivalent PNG and
SVG figures. The main figure covers sample medians, within-sample spread,
within-tissue similarity, and paired-sample correlation.

## 9. Automated validation

Twelve offline checks validate argument defaults, pinned source metadata, matrix
alignment, duplicate and missing-value detection, threshold calculation, sample
metrics, and participant-pair reconstruction.

## 10. Limitations

- Processed values depend on the original study's preprocessing decisions.
- Correlation flags do not identify the cause of unusual profiles.
- Retaining flagged arrays does not assert that all arrays are technically ideal.
- QC results do not establish differential expression or clinical significance.

## References

- [GEO series GSE19804](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE19804)
- [GPL570 platform](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL570)
- [NCBI GEO documentation](https://www.ncbi.nlm.nih.gov/geo/info/overview.html)
- [NCBI data-usage policies](https://www.ncbi.nlm.nih.gov/home/about/policies/)
