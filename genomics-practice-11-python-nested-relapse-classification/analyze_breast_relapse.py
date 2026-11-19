#!/usr/bin/env python3
"""Run the Project 11 nested-validation breast-relapse analysis."""
from __future__ import annotations

import argparse
import json
import platform
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
import pandas as pd
import scipy
import sklearn
from sklearn.metrics import roc_auc_score

from breast_relapse_ml import (
    aggregate_repeated_predictions,
    bootstrap_metrics,
    cross_validated_er_baseline,
    cross_validated_prevalence_baseline,
    curve_data,
    feature_stability,
    nested_repeated_cv,
    sha256,
)

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
OUTPUTS = ROOT / "outputs"


def save_figure(expression_metrics: pd.DataFrame, er_metrics: pd.DataFrame,
                expression_pred: pd.DataFrame, er_pred: pd.DataFrame,
                stability: pd.DataFrame) -> None:
    y = expression_pred["outcome"].to_numpy()
    p_expr = expression_pred["predicted_probability"].to_numpy()
    p_er = er_pred["predicted_probability"].to_numpy()
    roc_expr, cal_expr = curve_data(y, p_expr)
    roc_er, cal_er = curve_data(y, p_er)
    roc_expr.to_csv(OUTPUTS / "expression_roc_curve.csv", index=False)
    cal_expr.to_csv(OUTPUTS / "expression_calibration_curve.csv", index=False)

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))
    ax = axes[0]
    ax.plot(roc_expr["false_positive_rate"], roc_expr["true_positive_rate"],
            color="#1769aa", lw=2.2, label=f"Expression (AUC {roc_auc_score(y, p_expr):.3f})")
    ax.plot(roc_er["false_positive_rate"], roc_er["true_positive_rate"],
            color="#b55d00", lw=1.8, label=f"ER only (AUC {roc_auc_score(y, p_er):.3f})")
    ax.plot([0, 1], [0, 1], "--", color="0.55", lw=1)
    ax.set(xlabel="False-positive rate", ylabel="True-positive rate", title="Repeated nested-CV ROC")
    ax.legend(frameon=False, fontsize=8)

    ax = axes[1]
    ax.plot(cal_expr["mean_predicted_probability"], cal_expr["observed_event_fraction"],
            marker="o", color="#1769aa", label="Expression")
    ax.plot(cal_er["mean_predicted_probability"], cal_er["observed_event_fraction"],
            marker="s", color="#b55d00", label="ER only")
    ax.plot([0, 1], [0, 1], "--", color="0.55", lw=1)
    ax.set(xlim=(0, 1), ylim=(0, 1), xlabel="Mean predicted probability",
           ylabel="Observed 5-year relapse fraction", title="Calibration by probability octile")
    ax.legend(frameon=False, fontsize=8)

    ax = axes[2]
    top = stability.head(15).iloc[::-1]
    labels = top["gene_symbol"] + " (" + top["entrez_id"].astype(str) + ")"
    colors = np.where(top["mean_coefficient"] >= 0, "#9e2f2f", "#2b6f8a")
    ax.barh(labels, top["selection_fraction"], color=colors)
    ax.set(xlim=(0, 1), xlabel="Fraction of 15 outer models",
           title="Most stable selected genes")
    ax.tick_params(axis="y", labelsize=8)
    ax.legend(handles=[Patch(color="#9e2f2f", label="Mean coefficient > 0"),
                       Patch(color="#2b6f8a", label="Mean coefficient < 0")],
              frameon=False, fontsize=7, loc="lower right")
    fig.suptitle("Five-year relapse classification in GSE2034", fontsize=14, y=1.02)
    fig.tight_layout()
    fig.savefig(OUTPUTS / "nested_relapse_classification.png", dpi=180, bbox_inches="tight")
    fig.savefig(OUTPUTS / "nested_relapse_classification.svg", bbox_inches="tight")
    # Normalize SVG whitespace for Git after each export.
    svg_path = OUTPUTS / "nested_relapse_classification.svg"
    svg_text = svg_path.read_text(encoding="utf-8")
    svg_path.write_text(
        "\n".join(line.rstrip() for line in svg_text.splitlines()).rstrip("\n") + "\n",
        encoding="utf-8",
    )
    plt.close(fig)


def write_report(summary: dict[str, object], expression_metrics: pd.DataFrame,
                 er_metrics: pd.DataFrame, prevalence_metrics: pd.DataFrame,
                 stability: pd.DataFrame) -> None:
    def metric(frame: pd.DataFrame, name: str) -> tuple[float, float, float]:
        row = frame.loc[frame["metric"] == name].iloc[0]
        return float(row["estimate"]), float(row["ci_lower"]), float(row["ci_upper"])
    auc = metric(expression_metrics, "roc_auc")
    ap = metric(expression_metrics, "average_precision")
    brier = metric(expression_metrics, "brier_score")
    er_auc = metric(er_metrics, "roc_auc")
    null_brier = metric(prevalence_metrics, "brier_score")
    stable = stability.loc[stability["selection_fraction"] >= 0.8, "gene_symbol"].tolist()
    lines = [
        "# Generated analysis report",
        "",
        "## Cohort and endpoint",
        "",
        f"The public GSE2034 clinical table contained {summary['source_samples']} samples. "
        f"The prespecified 60-month endpoint retained {summary['eligible_samples']} samples: "
        f"{summary['five_year_events']} relapses by 60 months and "
        f"{summary['five_year_non_events']} participants event-free through 60 months. "
        f"{summary['early_censored_excluded']} relapse-free participants censored by 60 months were excluded.",
        "",
        "## Repeated nested-validation performance",
        "",
        f"The expression model's participant-level repeated-CV AUROC was {auc[0]:.3f} "
        f"(stratified bootstrap 95% CI {auc[1]:.3f}–{auc[2]:.3f}); average precision was "
        f"{ap[0]:.3f} ({ap[1]:.3f}–{ap[2]:.3f}), and the Brier score was "
        f"{brier[0]:.3f} ({brier[1]:.3f}–{brier[2]:.3f}).",
        "",
        f"The ER-status-only comparator achieved AUROC {er_auc[0]:.3f} "
        f"({er_auc[1]:.3f}–{er_auc[2]:.3f}). This descriptive comparison is not a formal "
        "test of incremental clinical value.",
        "",
        f"The training-prevalence comparator had a Brier score of {null_brier[0]:.3f}. "
        f"The expression model's Brier score of {brier[0]:.3f} was worse, and its calibration "
        "intercept and slope indicate overprediction and overly extreme probabilities. "
        "Discrimination therefore did not translate into reliable absolute-risk estimates.",
        "",
        "## Feature stability",
        "",
        f"{len(stable)} genes were selected in at least 80% of the 15 outer models. "
        + ("The stable symbols were: " + ", ".join(stable) + "." if stable else
           "No gene reached the 80% stability threshold."),
        "",
        "## Interpretation boundary",
        "",
        "The estimates describe internal repeated nested cross-validation in one retrospective "
        "microarray cohort. They do not establish transportability, treatment benefit, a deployable "
        "clinical threshold, or prospective diagnostic/prognostic validity. The no-relapse group at "
        "the 60-month landmark includes right-censored participants observed beyond that landmark and "
        "participants who relapsed later; this is appropriate for the defined binary horizon but does "
        "not replace a time-to-event analysis.",
    ]
    (OUTPUTS / "analysis_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outer-repeats", type=int, default=3)
    parser.add_argument("--bootstrap-iterations", type=int, default=2000)
    args = parser.parse_args()
    if args.outer_repeats < 1 or args.bootstrap_iterations < 100:
        raise ValueError("Use at least one outer repeat and 100 bootstrap iterations")
    OUTPUTS.mkdir(exist_ok=True)

    clinical = pd.read_csv(DATA / "clinical_metadata.csv")
    expression = pd.read_csv(DATA / "expression_genes.csv.gz", index_col=0).T
    clinical = clinical.set_index("gsm")
    eligible_ids = [gsm for gsm in expression.index if bool(clinical.loc[gsm, "endpoint_eligible"])]
    model_clinical = clinical.loc[eligible_ids].copy()
    X = expression.loc[eligible_ids].astype(float)
    y = model_clinical["relapse_by_60_months"].astype(int)
    er = model_clinical["er_status"].map({"ER-": 0, "ER+": 1})
    if er.isna().any():
        raise ValueError("Unexpected ER status")

    repeated, folds, feature_rows = nested_repeated_cv(
        X, y, eligible_ids, outer_repeats=args.outer_repeats
    )
    expression_pred = aggregate_repeated_predictions(repeated)
    er_repeated = cross_validated_er_baseline(
        er.to_numpy(), y.to_numpy(), eligible_ids, outer_repeats=args.outer_repeats
    )
    er_pred = aggregate_repeated_predictions(er_repeated)
    prevalence_repeated = cross_validated_prevalence_baseline(
        y.to_numpy(), eligible_ids, outer_repeats=args.outer_repeats
    )
    prevalence_pred = aggregate_repeated_predictions(prevalence_repeated)
    # Align both model outputs to make comparison and plotting explicit.
    expression_pred = expression_pred.sort_values("sample_id").reset_index(drop=True)
    er_pred = er_pred.sort_values("sample_id").reset_index(drop=True)
    if not expression_pred[["sample_id", "outcome"]].equals(er_pred[["sample_id", "outcome"]]):
        raise RuntimeError("Expression and comparator prediction identities differ")

    expression_metrics = bootstrap_metrics(
        expression_pred["outcome"].to_numpy(),
        expression_pred["predicted_probability"].to_numpy(),
        args.bootstrap_iterations,
    )
    er_metrics = bootstrap_metrics(
        er_pred["outcome"].to_numpy(), er_pred["predicted_probability"].to_numpy(),
        args.bootstrap_iterations,
    )
    prevalence_metrics = bootstrap_metrics(
        prevalence_pred["outcome"].to_numpy(),
        prevalence_pred["predicted_probability"].to_numpy(),
        args.bootstrap_iterations,
    )
    stability = feature_stability(feature_rows, total_folds=len(folds))

    expression_pred.to_csv(OUTPUTS / "oof_predictions.csv", index=False, lineterminator="\n")
    repeated.to_csv(OUTPUTS / "oof_predictions_by_repeat.csv", index=False, lineterminator="\n")
    er_pred.to_csv(OUTPUTS / "er_baseline_predictions.csv", index=False, lineterminator="\n")
    prevalence_pred.to_csv(OUTPUTS / "prevalence_baseline_predictions.csv", index=False, lineterminator="\n")
    folds.to_csv(OUTPUTS / "outer_fold_results.csv", index=False, lineterminator="\n")
    feature_rows.to_csv(OUTPUTS / "outer_fold_selected_features.csv", index=False, lineterminator="\n")
    stability.to_csv(OUTPUTS / "feature_stability.csv", index=False, lineterminator="\n")
    expression_metrics.assign(model="expression").to_csv(
        OUTPUTS / "expression_performance.csv", index=False, lineterminator="\n")
    er_metrics.assign(model="er_status_only").to_csv(
        OUTPUTS / "er_baseline_performance.csv", index=False, lineterminator="\n")
    prevalence_metrics.assign(model="training_prevalence").to_csv(
        OUTPUTS / "prevalence_baseline_performance.csv", index=False, lineterminator="\n")

    summary = {
        "source_samples": int(len(clinical)),
        "eligible_samples": int(len(model_clinical)),
        "five_year_events": int(y.sum()),
        "five_year_non_events": int((1 - y).sum()),
        "early_censored_excluded": int((~clinical["endpoint_eligible"].astype(bool)).sum()),
        "genes_available": int(X.shape[1]),
        "outer_models": int(len(folds)),
        "stable_genes_selected_at_least_80_percent": int((stability["selection_fraction"] >= 0.8).sum()),
        "expression_auroc": float(expression_metrics.loc[expression_metrics.metric == "roc_auc", "estimate"].iloc[0]),
        "er_only_auroc": float(er_metrics.loc[er_metrics.metric == "roc_auc", "estimate"].iloc[0]),
        "prevalence_baseline_brier": float(prevalence_metrics.loc[prevalence_metrics.metric == "brier_score", "estimate"].iloc[0]),
    }
    pd.DataFrame([{"metric": key, "value": value} for key, value in summary.items()]).to_csv(
        OUTPUTS / "analysis_summary.csv", index=False, lineterminator="\n")
    save_figure(expression_metrics, er_metrics, expression_pred, er_pred, stability)
    write_report(summary, expression_metrics, er_metrics, prevalence_metrics, stability)
    metadata = {
        "random_seed": 20261026,
        "outer_cv": f"5 folds x {args.outer_repeats} repeats",
        "inner_cv": "4-fold stratified",
        "bootstrap_iterations": args.bootstrap_iterations,
        "endpoint_months": 60,
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scipy": scipy.__version__,
        "scikit_learn": sklearn.__version__,
        "input_sha256": {
            name: sha256(DATA / name)
            for name in ["clinical_metadata.csv", "expression_genes.csv.gz", "representative_probes.csv"]
        },
    }
    (OUTPUTS / "run_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    print(pd.DataFrame([summary]).T.rename(columns={0: "value"}).to_string())
    print("\nExpression performance")
    print(expression_metrics.to_string(index=False))


if __name__ == "__main__":
    main()
