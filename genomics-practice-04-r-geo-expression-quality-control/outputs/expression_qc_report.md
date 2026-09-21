# Paired GEO Expression-Matrix Quality Assessment

## Dataset

- GEO series: `GSE19804`
- Platform: `GPL570`
- Design: paired lung tumor and adjacent-normal tissues
- Processed values: authors' quantile-normalized log2 expression matrix

## Results

- Probe sets: 54,675
- Samples: 120
- Complete tumor-normal pairs: 60
- Expression range: 3.0427 to 14.8903
- Sample-median range: 6.02783 to 6.02919
- Median within-pair correlation: 0.9489
- Samples marked for review: 8
- Samples automatically excluded: 0

## Samples marked for review

| GEO accession | Sample | Tissue | Mean within-tissue correlation |
| --- | --- | --- | ---: |
| GSM494571 | Lung Cancer 103T | Tumor | 0.8815 |
| GSM494572 | Lung Cancer 106T | Tumor | 0.8865 |
| GSM494582 | Lung Cancer 122T | Tumor | 0.8946 |
| GSM494596 | Lung Cancer 137T | Tumor | 0.8403 |
| GSM494654 | Lung Normal 135N | Normal | 0.8722 |
| GSM494656 | Lung Normal 137N | Normal | 0.9246 |
| GSM494657 | Lung Normal 139N | Normal | 0.8704 |
| GSM494675 | Lung Normal 189N | Normal | 0.9094 |

## Interpretation

The narrow sample-median range is consistent with the normalization reported
for the processed GEO matrix. Review flags identify unusually low mean
correlation within a tissue group using a Q1 minus 2 x IQR rule. Flagged
samples remain in the dataset because the differences may reflect technical
variation, participant-specific biology, or tumor heterogeneity.

Quality-control patterns do not establish differential expression, disease
causation, diagnostic performance, or clinical validity.
