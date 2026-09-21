# Analysis Workflow: Paired Lung-Tissue Differential Expression

## 1. Objective

Identify reproducible expression differences between matched lung-tumor and
adjacent-normal tissues while preserving participant pairing and documenting
sensitivity to sample-level QC flags.

## 2. Source verification

The complete `GSE19804` processed expression matrix and official `GPL570`
annotation are accepted only when their pinned byte counts and MD5 checksums
match. Source files are cached under `data/` and excluded from Git.

## 3. Pair reconstruction

GEO titles encode a participant number and tissue suffix. The workflow verifies
exactly one tumor and one adjacent-normal sample for each of 60 participants and
orders both matrices by participant before testing.

## 4. Primary analysis

For every probe set:

1. Calculate tumor-minus-normal expression within each participant.
2. Estimate the mean paired log2 fold change.
3. Perform a two-sided paired t-test across 60 differences.
4. Adjust 54,675 p-values using the Benjamini–Hochberg procedure.
5. Require FDR ≤ 0.05 and absolute mean log2 fold change ≥ 1.

## 5. Sensitivity analysis

Project 04's correlation rule is reconstructed directly from the matrix. Eight
arrays are marked for review across seven participants. Those seven complete
pairs are omitted, leaving 53 pairs.

A robust probe must meet the same FDR and effect-size criteria, with the same
direction, in both primary and sensitivity analyses.

## 6. Annotation integration

Probe identifiers are joined to GPL570 gene symbols, Entrez Gene identifiers,
and gene titles. Multi-gene symbols separated by `///` are retained in the full
and probe-level results but excluded from the unique unambiguous-gene summary.

When multiple robust probes map unambiguously to the same gene, the unique-gene
table retains the probe with the largest absolute primary effect, breaking ties
with the smaller FDR.

## 7. Outputs

The analysis writes the complete probe-level statistics, robust probe results,
unique-gene results, pairing and QC manifests, a report, and PNG/SVG figures.

## 8. Automated validation

Eleven offline tests cover pinned source metadata, MD5 calculation,
Benjamini–Hochberg adjustment, participant matching, sensitivity exclusion,
paired-test shape and direction, symbol ambiguity, and dimension validation.

## 9. Limitations

- The paired t-test assumes participant-level differences are adequately
  summarized by their mean and standard deviation.
- The selected effect-size and FDR thresholds are analytical criteria, not
  clinical decision boundaries.
- GPL570 annotation is platform-release specific and contains ambiguous mappings.
- The sensitivity analysis evaluates robustness to one declared QC rule; it does
  not exhaust all possible preprocessing decisions.

## References

- [GEO series GSE19804](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE19804)
- [GPL570 platform](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL570)
- [Originating study: Lu et al.](https://pubmed.ncbi.nlm.nih.gov/20802022/)
- [NCBI data-usage policies](https://www.ncbi.nlm.nih.gov/home/about/policies/)
