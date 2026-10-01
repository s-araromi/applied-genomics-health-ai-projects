# CFTR ClinVar Classification-Conflict Audit

**Snapshot date:** 2026-09-28

## Scope

The exact ClinVar query returned **549** records; **407** had the exact aggregate germline label `Conflicting classifications of pathogenicity` in ESummary. The audit examined the **25** eligible records with the largest numbers of supporting SCV accessions. Every selected VCV archive resolved as `current` during retrieval.

## Verified results

- Current records: **25/25**
- Records with a current GRCh38 coordinate: **25/25**
- Polarized pathogenic-tier versus benign-tier disagreements: **7**
- Records with criteria provided: **25**
- Records last evaluated more than five years before the snapshot: **0**
- Median aggregate classification count: **25**
- Maximum aggregate classification count: **35**

## Priority tiers

- Tier 1: polarized, multi-submitter: **7**
- Tier 2: criteria-provided, multi-submitter: **18**

Tier 1 requires a pathogenic-tier/benign-tier polarization, criteria-provided review status, and at least three submitters. Tier 2 requires criteria-provided review status and at least three submitters without that polarization. Tier 3 contains records with more limited review support. These are workflow priorities for manual evidence reconciliation, not variant classifications.

## Highest-submission records

- [VCV000007229.102](https://www.ncbi.nlm.nih.gov/clinvar/variation/7229/): 35 submissions from 33 submitters; Tier 1: polarized, multi-submitter.
- [VCV000053460.80](https://www.ncbi.nlm.nih.gov/clinvar/variation/53460/): 34 submissions from 31 submitters; Tier 2: criteria-provided, multi-submitter.
- [VCV000035835.90](https://www.ncbi.nlm.nih.gov/clinvar/variation/35835/): 31 submissions from 29 submitters; Tier 1: polarized, multi-submitter.
- [VCV000219537.98](https://www.ncbi.nlm.nih.gov/clinvar/variation/219537/): 31 submissions from 30 submitters; Tier 2: criteria-provided, multi-submitter.
- [VCV000007165.98](https://www.ncbi.nlm.nih.gov/clinvar/variation/7165/): 30 submissions from 27 submitters; Tier 1: polarized, multi-submitter.

## Interpretation

A large submission count does not make one interpretation correct. It indicates a larger evidence-reconciliation workload. Polarized records deserve particular attention because the aggregate conflict crosses the pathogenic-to-benign boundary. Clinical interpretation still requires condition-specific review of individual SCV assertions, inheritance, phenotype, allele frequency, functional evidence, and the most recent expert guidance.

## Limitations

- This is a high-submission subset, not all conflicting CFTR records.
- ClinVar is dynamic; later reruns can change record versions and rankings.
- Aggregate counts can combine condition-specific assertions and should not be treated as votes.
- The priority tiers are transparent operational rules, not validated clinical decision rules.
- No patient-level data were used, and the outputs must not be used for diagnosis or treatment.

## Source

- [ClinVar query](https://www.ncbi.nlm.nih.gov/clinvar/?term=CFTR[gene]+AND+"conflicting+classifications+of+pathogenicity"[clinsig])
- [ClinVar programmatic access](https://www.ncbi.nlm.nih.gov/clinvar/docs/programmatic_access/)
- [ClinVar review status](https://www.ncbi.nlm.nih.gov/clinvar/docs/review_status/)
