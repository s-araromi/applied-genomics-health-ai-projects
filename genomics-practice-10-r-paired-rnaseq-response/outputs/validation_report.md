# Validation Record

Verification completed on **19 October 2026** for the delivered Project 10.
Results below describe executed checks, not expected outcomes.

| Check | Executed result |
| --- | --- |
| Live source assessment | GSE52778 public; correct accession, study type, and all 16 expected series sample IDs |
| Count archive | airway 1.32.0 / Bioconductor 3.23; 24,756,071 bytes; pinned SHA-256 verified |
| Source reconstruction | Fresh retrieval and export reproduced `counts.csv`, `samples.csv`, and `gene_annotation.csv` hashes |
| Input coverage | 63,677 unique gene IDs; eight libraries; four complete donor pairs; no missing counts |
| R automated tests | 23 passed; complete results in `test_results.txt` |
| Model fitting | One paired model, one unpaired model, four donor-omission models; all parametric dispersion fits and coefficients converged |
| Rerun agreement | All 17 numerical CSV files reproduced byte-for-byte in a separate output directory |
| Independent numerical audit | Python/SciPy reproduced median-of-ratios factors, both Wald p-value families, threshold-test BH q values, and selection counts |
| Figures | Both PNG figures inspected; matching SVG versions included |
| R syntax | Every delivered R source parsed successfully |
| Safe import | Read-only default, missing Project 09 refusal, pending-change refusal, incorrect-origin refusal, modified-package refusal, successful apply, duplicate-project refusal |
| Git integrity | Existing project files preserved; staged whitespace check passed; committed-source hashes survived Git text normalization |

## Independent numerical agreement

The optional audit in `tests/audit_numerics.py` uses NumPy and SciPy rather than
DESeq2's result or normalization implementations. It independently computes
geometric means, median log ratios, normal survival probabilities, and BH
adjustment. Maximum absolute discrepancies were:

| Quantity | Maximum difference |
| --- | ---: |
| Size factors | 3.77e-15 |
| Point-null p values | 5.00e-15 |
| Effect-threshold p values | 5.50e-15 |
| Effect-threshold BH q values | 1.44e-15 |

Relative tolerances also covered very small nonzero probabilities; the audit
used Python's full floating-point parsing. This validates the numerical
translation of the chosen model outputs. It does not independently establish
that the negative-binomial model is the true biological data-generating model.

## Verified findings

The predeclared count filter retained 16,139 genes. The full paired model
selected 600 genes at effect-threshold FDR < 0.05: 352 increases and 248
decreases. All 600 had the same effect direction in all four donor-omission
refits. Of these, 359 remained significant in all four. The point-null family
selected 4,068 genes. The unpaired effect-threshold model selected 292 genes,
with 289 shared discoveries. These summaries were reconciled to the complete
exported tables.

## Boundaries

The verified software environment was R 4.3.3, DESeq2 1.42.0,
SummarizedExperiment 1.32.0, and digest 0.6.34 on Ubuntu 24.04.3. All
loaded dependency versions are in `session_info.txt`. The data distribution
was airway 1.32.0; it was loaded directly, not installed as an analysis
software dependency. The install script selects packages compatible with the
user's R version; a newer DESeq2 release has not been claimed to reproduce
these results exactly.

Importer checks used an isolated local structural repository with the actual
Projects 01–08 and minimal Project 09 folder metadata, because the delivered
Project 09 was not yet on remote main. This checks file and Git behaviour;
it does not represent a new validation of Project 09 or write to the user's
repository. No commit or push to GitHub was made.

Raw-read QC, alignment, re-counting, independent biological replication,
predictive performance, and clinical efficacy were outside the performed
work. The overlapping donor-omission fits are sensitivity analyses. No
simulated or substituted biological dataset was used.
