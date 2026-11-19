#!/usr/bin/env python3
"""Independent consistency checks for generated Project 11 outputs."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "outputs"


def summary_value(summary: pd.DataFrame, name: str) -> float:
    return float(summary.loc[summary["metric"] == name, "value"].iloc[0])


def main() -> None:
    summary = pd.read_csv(OUT / "analysis_summary.csv")
    pred = pd.read_csv(OUT / "oof_predictions.csv")
    repeated = pd.read_csv(OUT / "oof_predictions_by_repeat.csv")
    folds = pd.read_csv(OUT / "outer_fold_results.csv")
    selected = pd.read_csv(OUT / "outer_fold_selected_features.csv")
    stability = pd.read_csv(OUT / "feature_stability.csv")
    performance = pd.read_csv(OUT / "expression_performance.csv")
    metadata = json.loads((OUT / "run_metadata.json").read_text())

    assert len(pred) == 276 and pred["sample_id"].is_unique
    assert pred["outcome"].sum() == 95
    assert len(repeated) == 828
    assert (repeated.groupby("sample_id").size() == 3).all()
    assert len(folds) == 15 and set(folds["best_k"]).issubset({25, 50, 100})
    selected_counts = selected.groupby(["repeat", "fold"]).size().sort_index()
    fold_counts = folds.set_index(["repeat", "fold"])["best_k"].sort_index()
    assert selected_counts.equals(fold_counts)
    assert (stability["selection_fraction"] <= 1).all()
    assert (stability["selection_count"] <= 15).all()

    y = pred["outcome"].to_numpy()
    p = pred["predicted_probability"].to_numpy()
    recomputed = {
        "roc_auc": roc_auc_score(y, p),
        "average_precision": average_precision_score(y, p),
        "brier_score": brier_score_loss(y, p),
    }
    for name, value in recomputed.items():
        recorded = float(performance.loc[performance["metric"] == name, "estimate"].iloc[0])
        assert np.isclose(value, recorded, rtol=0, atol=1e-12)
    assert np.isclose(recomputed["roc_auc"], summary_value(summary, "expression_auroc"), atol=1e-12)
    assert summary_value(summary, "eligible_samples") == 276
    assert summary_value(summary, "stable_genes_selected_at_least_80_percent") == int(
        (stability["selection_fraction"] >= 0.8).sum()
    )
    assert metadata["endpoint_months"] == 60
    assert metadata["outer_cv"] == "5 folds x 3 repeats"
    print("Verified prediction coverage, nested-fold accounting, metrics, stability, and run metadata.")


if __name__ == "__main__":
    main()
