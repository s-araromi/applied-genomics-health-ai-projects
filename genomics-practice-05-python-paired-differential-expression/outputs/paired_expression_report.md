# Paired Differential-Expression Analysis

## Dataset

- GEO series: `GSE19804`
- Platform: `GPL570`
- Primary analysis: 60 paired tumor-normal participants
- Sensitivity analysis: 53 pairs

## Results

- Probe sets tested: 54,675
- Primary signals meeting FDR and effect-size criteria: 1,960
- Higher in tumor: 591
- Lower in tumor: 1,369
- Robust probe sets after sensitivity analysis: 1,928
- Robust unique, unambiguous genes: 1,328
- Probe sets with gene symbols: 45,118
- Samples marked for review: 8

## Leading robust signals

Higher in tumor: GOLM1, OCIAD2, PYCR1, COL10A1, ALDH18A1

Lower in tumor: SPTBN1, AGER, CA4, SPOCK2, GPM6A

Three SEMA5A probes were robustly lower in tumor, with primary mean
log2 fold changes from -1.44 to -1.78. This agrees with the originating
study's interest in SEMA5A but is not an independent replication because
the same cohort is analyzed here.

## Interpretation

Signals were required to meet BH FDR <= 0.05 and an absolute paired
mean log2 fold change >= 1. Robust signals also met both criteria with
the same direction after excluding seven pairs containing review-flagged
arrays. Association in this cohort does not establish causation, diagnostic
performance, or clinical utility.
