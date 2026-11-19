# Generated analysis report

## Cohort and endpoint

The public GSE2034 clinical table contained 286 samples. The prespecified 60-month endpoint retained 276 samples: 95 relapses by 60 months and 181 participants event-free through 60 months. 10 relapse-free participants censored by 60 months were excluded.

## Repeated nested-validation performance

The expression model's participant-level repeated-CV AUROC was 0.650 (stratified bootstrap 95% CI 0.581–0.713); average precision was 0.487 (0.415–0.572), and the Brier score was 0.235 (0.215–0.256).

The ER-status-only comparator achieved AUROC 0.476 (0.404–0.547). This descriptive comparison is not a formal test of incremental clinical value.

The training-prevalence comparator had a Brier score of 0.226. The expression model's Brier score of 0.235 was worse, and its calibration intercept and slope indicate overprediction and overly extreme probabilities. Discrimination therefore did not translate into reliable absolute-risk estimates.

## Feature stability

11 genes were selected in at least 80% of the 15 outer models. The stable symbols were: NEK2, RACGAP1, EFR3A, TTI1, CCNE2, TCEB1, SNU13, CD74, FRAS1, GINS1, KPNA2.

## Interpretation boundary

The estimates describe internal repeated nested cross-validation in one retrospective microarray cohort. They do not establish transportability, treatment benefit, a deployable clinical threshold, or prospective diagnostic/prognostic validity. The no-relapse group at the 60-month landmark includes right-censored participants observed beyond that landmark and participants who relapsed later; this is appropriate for the defined binary horizon but does not replace a time-to-event analysis.
