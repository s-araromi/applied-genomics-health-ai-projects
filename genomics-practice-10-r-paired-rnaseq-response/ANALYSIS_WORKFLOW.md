# Project 10 Analytical Workflow and Reproducibility Record

## Scope and analytical unit

The analysis concerns gene-level relative fragment abundance in eight RNA-seq
libraries from four donor-derived airway smooth-muscle cell lines. A biological
replicate is a donor, not a gene, read, or pair of sequencing mates. The target
contrast is dexamethasone versus untreated cells from the same donor. Only the
eight libraries with no albuterol treatment were included. All 63,677 gene
rows in the selected source object were retained during export.

The work began from an authentic, previously counted matrix. It did not perform
FASTQ quality control, alignment, fragment counting, or an audit of read-level
identity. The original study's FPKM supplementary table and the data package's
small example BAM subsets were not used in statistical modelling.

## Source verification and provenance

The source is airway 1.32.0 from Bioconductor release 3.23, commit deb0a4f.
`data/archive_pin.json` records the versioned URLs, archive length, SHA-256,
release, and original reference. The package includes the original public
GSE52778 RNA-seq experiment and credits Himes et al. The count object has
63,677 rows and these columns, in this order:

| Donor | Untreated run / GEO sample | Dexamethasone run / GEO sample |
| --- | --- | --- |
| N61311 | [SRR1039508](https://www.ncbi.nlm.nih.gov/sra?term=SRR1039508) / [GSM1275862](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1275862) | [SRR1039509](https://www.ncbi.nlm.nih.gov/sra?term=SRR1039509) / [GSM1275863](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1275863) |
| N052611 | [SRR1039512](https://www.ncbi.nlm.nih.gov/sra?term=SRR1039512) / [GSM1275866](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1275866) | [SRR1039513](https://www.ncbi.nlm.nih.gov/sra?term=SRR1039513) / [GSM1275867](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1275867) |
| N080611 | [SRR1039516](https://www.ncbi.nlm.nih.gov/sra?term=SRR1039516) / [GSM1275870](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1275870) | [SRR1039517](https://www.ncbi.nlm.nih.gov/sra?term=SRR1039517) / [GSM1275871](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1275871) |
| N061011 | [SRR1039520](https://www.ncbi.nlm.nih.gov/sra?term=SRR1039520) / [GSM1275874](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1275874) | [SRR1039521](https://www.ncbi.nlm.nih.gov/sra?term=SRR1039521) / [GSM1275875](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1275875) |

The live GEO series metadata were retrieved on 19 October 2026 using:

```text
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52778&targ=self&view=full&form=text
```

The response identified the series as public, listed the expected 16 original
samples, and showed a last-update date of 13 January 2023. The response, normalized to LF line endings,
is committed in `data/GSE52778_status_2026-10-19.soft`; it was not rewritten as
an inferred status. Archive availability and GEO status were assessed
separately. The initial primary archive request returned HTTP 403; retrieval
succeeded through the officially listed Posit mirror. A later end-to-end
retrieval test also succeeded against the primary URL with the same hash.

The count-reference provenance is GRCh37.p13 and Ensembl release 75. The airway
vignette documents STAR alignment and `summarizeOverlaps` using union mode,
paired fragments, and ignored strand. January 2023 gene labels in the object
were preserved. The old annotation is part of the source identity, not a claim
that these coordinates or symbols are current clinical references.

## Data flow and file responsibilities

1. `retrieve_data.R` downloads the exact archive and checks its SHA-256 before
   opening it. It tries only the primary endpoint and its official mirror;
   both must supply identical bytes. It extracts only the count object and
   package DESCRIPTION to temporary storage.
2. `R/export_airway.R` validates the selected object and exports `counts.csv`,
   `samples.csv`, and `gene_annotation.csv`. Every gene row is exported. The
   committed table hashes must be reproduced; a different package cannot
   silently replace them. The download destination must be empty.
3. `data/source_manifest.csv` pins those three tables, the package DESCRIPTION,
   and the live GEO status response. `read_inputs()` checks all five identities
   before fitting. The delivery manifest separately covers the complete package.
4. `R/rnaseq_functions.R` validates, aligns, filters, fits models, performs BH
   correction, and computes descriptive sample summaries. Sample metadata are
   matched by run ID; positional metadata assignment is prohibited.
5. `analyze_rnaseq.R` runs all six fits, exports results, and calls `R/plots.R`.
6. `tests/test_rnaseq.R` tests source identities, valid pairing, malformed-input
   rejection, testing-family reconciliation, PCA geometry, and exported results.
7. `verify_outputs.R` compares all numerical CSV files from a clean run against
   the delivered reference. Package versions and plots are reviewed separately.

## Input validation

Counts must be a nonempty matrix of finite, non-negative integers within the
R integer range. Gene and library IDs must be unique and nonmissing. All
libraries must have matching metadata; every donor must contribute exactly
one untreated and one dexamethasone sample, with no albuterol. At least three
complete donors are required. The design matrix must be full rank with positive
residual degrees of freedom. Annotation rows must align exactly to gene IDs.
No sample was discarded after observing analysis outcomes.

Assigned gene-fragment totals were 20,637,971; 18,809,481; 25,348,649;
15,163,415; 24,448,408; 30,818,215; 19,126,151; and 21,164,133 in source order.
They were checked against the pinned airway object. These totals are not
alignment counts for the entire libraries.

## Gene filtering and normalization

The prespecified filter requires count >= 10 in at least four libraries.
It retained 16,139 genes and excluded 47,538. Excluded genes are recorded in
`gene_filter.csv`; they were not assigned differential-expression p values.
The same retained-gene family was used in the unpaired and donor-omission
analyses, so changes in discovery counts do not arise from a new gene filter.
The filter uses all eight libraries and is part of a descriptive experiment
analysis, not a predictive-training procedure.

DESeq2 uses median-of-ratios normalization on the retained matrix. For genes
with positive counts in every library, compute the geometric mean across
libraries; each sample's size factor is the exponential of the median log
ratio to those means. `sample_qc.csv` reports the resulting factors and both
full-source and retained-gene fragment totals. Normalized counts equal counts
divided by size factors and still include donor differences.

Normalization assumes a sufficiently stable reference distribution. It does
not identify absolute global RNA shifts. No external spike-in calibration or
transcript-length correction was introduced. The analysis concerns relative
fragment abundance for the same gene across conditions, not comparisons of
expression magnitude between genes.

## Negative-binomial model

For gene i and library j:

```text
K_ij ~ NegativeBinomial(mu_ij, alpha_i)
mu_ij = s_j * q_ij
log2(q_ij) = intercept_i + donor_effect_i,j + beta_i * treated_j
Var(K_ij) = mu_ij + alpha_i * mu_ij^2
```

`~ donor + dex` represents donor through fixed intercept differences. There
are eight observations and five design columns, giving three residual degrees
of freedom. Untreated is explicitly the reference, so positive beta means
higher relative expression after dexamethasone. A separate interaction for
every donor-condition cell would saturate the model and is not fitted.

DESeq2 estimates gene dispersions, a parametric mean-dispersion trend, and
moderated gene-wise dispersions. `betaPrior=FALSE` keeps coefficients unshrunk;
no claim of shrunken effect-size ranking is made. `parallel=FALSE` avoids
worker-order variability. `minReplicatesForReplace=Inf` disables count
replacement. The code stops if coefficient fitting does not converge or if
the declared parametric dispersion fit changes. All six fits used the
parametric trend and converged in the verified execution.

`dispersion_estimates.csv` exposes gene estimates, trend values, and final
dispersions. `influence_diagnostics.csv` reports each gene's maximum Cook's
distance and the nominal 99th-percentile F(5,3) reference cutoff for the full
model. DESeq2's default Cook's filtering rules remain enabled in `results()`;
the exported cutoff is a diagnostic reference, not a manual deletion rule.
No primary p values were unavailable in this run. No sample or gene was
removed in response to an influence statistic.

## Hypotheses, intervals, and multiplicity

<!-- Project 10 threshold-method compatibility correction -->
### Threshold-test compatibility

Project 10 preserves the conservative 2014 threshold Wald test used in its
original DESeq2 1.42.0 analysis. DESeq2 1.44.0 changed the implementation behind
`greaterAbs`; this project explicitly uses `greaterAbs2014` on DESeq2 >= 1.44.0
and `greaterAbs` on earlier releases. The choice preserves the original
statistical procedure; it does not tune the analysis to a target gene count.
Other differences across software environments may still affect numerical values.

An independent calculation using the original fitted coefficients reproduces
600 primary genes and 359 selected in every donor-omission refit with the
original method, versus 687 and 396 with the newer method. The newer method
also yields 391 increases, 296 decreases, and 340 unpaired discoveries, matching
the reported DESeq2 1.52.0 rerun. These are different implementations of the
same interval-null hypothesis, rather than a change to the source data.

The test suite now includes an additional check that exported threshold
p-values match the original normal-tail formula across all six fitted models.
Run the analysis before the tests after applying this correction; the suite
should then report 24 tests. Original 23-test logs describe the original release.

Source: [DESeq2 release notes, version 1.44.0](https://github.com/thelovelab/DESeq2/blob/devel/NEWS).


The primary hypothesis is an interval null:

```text
H0: |beta_i| <= log2(1.5)
H1: |beta_i| > log2(1.5)
```

It was implemented with `lfcThreshold=log2(1.5)` and
`altHypothesis="greaterAbs2014"` on DESeq2 >= 1.44.0 (the original
`greaterAbs` implementation on earlier versions). This tests evidence that the population model
coefficient exceeds the effect threshold. It is not a point-null p value
followed by an observed-fold-change cutoff.

A secondary two-sided point-null test evaluates beta_i = 0. The two p-value
families are corrected separately and labelled explicitly. Independent
filtering is disabled in both. BH adjustment uses all 16,139 retained genes;
an unavailable p value is treated as one for adjustment and then its q value
is set back to NA. Missing tests cannot become discoveries. This conservative
fixed-family convention is also used within each sensitivity refit.

The primary selected set is exactly nonmissing `q_above_1_5_fold < 0.05`.
Fitted beta values and standard errors also define marginal 95% Wald intervals.
The intervals are not simultaneous and do not replace the primary threshold
test. FDR is a model-based expected proportion across a discovery set, not a
probability that any named gene is false.

## Dimensionality reduction and visualization

A design-aware variance-stabilizing transformation (`blind=FALSE`) reduces the
count mean-variance dependence for sample exploration. The 500 retained genes
with highest transformed variance are selected with gene-ID tie breaking.
`prcomp` operates on centered sample-by-gene values without scaling every gene
to unit variance. PC signs are oriented deterministically using the largest
absolute loading, because an eigenvector's sign is otherwise arbitrary.

PC1 and PC2 explain 48.576% and 23.274% of variance in that selected feature
matrix. Donor pairs are joined in the plot, with treatment encoded by color
and donor by shape. This is a descriptive projection of the same eight samples,
using full-data normalization and design-aware transformation; separation is
not predictive accuracy or independent validation.

`sample_distances.csv` uses Euclidean distance across all retained VST genes.
The top-gene heatmap uses log2[(normalized treated + 0.5)/(normalized untreated
+ 0.5)] for each donor. It displays exact numbers and a common color scale;
these ratios are descriptive, pseudocount dependent, and distinct from the
negative-binomial coefficients.

## Design and donor sensitivity

The unpaired model `~ dex` estimates the same treatment contrast while treating
eight libraries as independent. It has six residual degrees of freedom but
ignores the matched biological structure. It is a sensitivity comparator,
not the preferred analysis or a competing model chosen for the largest list.
It selected 292 genes, of which 289 overlapped the paired selection. Thus 311
paired discoveries were not selected by the unpaired model and three unpaired
discoveries were not selected by the paired model.

Each leave-one-donor-out fit removes both libraries from one donor, refits
size factors and dispersions, and uses `~ donor + dex` on six libraries from
three donors. These models have two residual degrees of freedom. The same
16,139 gene hypotheses remain in each refit. Per-gene exports include fitted
effects, threshold tests, number of matching signs, number of significant
refits, smallest absolute omission effect, and largest absolute effect change.

All 600 full-model primary genes retained their direction in all four refits.
Counts significant in zero through four refits were 34, 98, 71, 38, and 359.
The four refits overlap heavily; no independent replication rate or new FDR
claim is assigned to this stability summary. Significance loss can arise from
reduced sample size as well as influence. A significant full-model finding
was not retroactively redefined by the sensitivity results.

## Reproduction commands

Run from this project directory. No Python runtime is required for the main
analysis. R 4.3+ and compatible DESeq2/SummarizedExperiment/digest packages are
required. Network access is needed for package installation and optional source
retrieval, but not for the bundled analysis.

```bash
Rscript --vanilla install_dependencies.R
Rscript --vanilla tests/test_rnaseq.R
Rscript --vanilla analyze_rnaseq.R
```

To verify source export independently, choose an empty destination:

```bash
Rscript --vanilla retrieve_data.R retrieved_data
```

To reproduce without replacing the delivered outputs:

```bash
Rscript --vanilla analyze_rnaseq.R rerun_outputs
Rscript --vanilla verify_outputs.R outputs rerun_outputs
```

Add `--strict` to the comparison command for identical CSV bytes. The ordinary
comparison checks table structure, identifiers, and numeric values with
`all.equal` tolerance 1e-7. Byte identity was verified in the original runtime;
newer DESeq2, R, BLAS, or platform versions may change values. A tolerance
failure requires investigation rather than automatic output replacement.
`outputs/session_info.txt` records the verified environment. Installation does
not claim to recreate an immutable software container.

## Validation and completed checks

The delivered execution passed 23 tests. A clean rerun reproduced every
numerical CSV byte-for-byte. Source retrieval reproduced the three input-table
SHA-256 values from a fresh archive download. An independent SciPy calculation
reconciled the R size factors, point-null and effect-threshold Wald p values,
and BH-adjusted q values. Both figures were visually inspected. An isolated
Git import was checked for pending-change refusal, source-package corruption,
duplicate project protection, preservation of existing README content, and
clean whitespace. See `outputs/validation_report.md` for exact outcomes.

## Scientific limitations and usage boundaries

The donor sample is small and demographically narrow. One dose and time point
in primary cell culture cannot establish clinical benefit or transportability
to patients. The fitted treatment effect is shared across donors and does not
estimate population heterogeneity reliably. Count processing uses an older
reference and may differ from modern quantification. Relative normalization
cannot recover absolute transcript changes without external calibration.
Results depend on model assumptions and threshold choice. No predictive model,
external cohort validation, read-level QC, new patient dataset, or mechanistic
causal analysis was performed. Public source terms, author attribution, and
non-identification boundaries are recorded in `DATA_USAGE.md`.
