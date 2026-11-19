"""Core data and modelling functions for Project 11.

The module keeps source parsing, endpoint construction, nested validation, and
metric calculation testable independently from the command-line workflow.
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import io
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from scipy.special import expit, logit
from sklearn.calibration import calibration_curve
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    average_precision_score,
    balanced_accuracy_score,
    brier_score_loss,
    confusion_matrix,
    roc_auc_score,
    roc_curve,
)
from sklearn.model_selection import GridSearchCV, RepeatedStratifiedKFold, StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


RANDOM_SEED = 20261026


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _split_geo_line(line: str) -> list[str]:
    return next(csv.reader([line], delimiter="\t", quotechar='"'))


def parse_series_matrix(path: Path) -> tuple[pd.DataFrame, dict[str, list[str]]]:
    """Return probe-by-sample values and selected GEO sample metadata."""
    metadata: dict[str, list[str]] = {}
    table_lines: list[str] = []
    in_table = False
    with gzip.open(path, "rt", encoding="utf-8", errors="strict") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line.startswith("!Sample_"):
                fields = _split_geo_line(line)
                metadata.setdefault(fields[0], fields[1:])
            elif line == "!series_matrix_table_begin":
                in_table = True
            elif line == "!series_matrix_table_end":
                break
            elif in_table:
                table_lines.append(line)
    if not table_lines:
        raise ValueError("GEO Series Matrix expression table was not found")
    frame = pd.read_csv(io.StringIO("\n".join(table_lines)), sep="\t", index_col=0)
    frame.index = frame.index.astype(str)
    frame.columns = frame.columns.astype(str)
    if frame.index.has_duplicates or frame.columns.has_duplicates:
        raise ValueError("Probe or sample identifiers are not unique")
    if not np.isfinite(frame.to_numpy(dtype=float)).all():
        raise ValueError("Expression matrix contains non-finite values")
    return frame, metadata


def parse_clinical_table(path: Path) -> pd.DataFrame:
    """Extract the patient clinical table from the GSE2034 family SOFT file."""
    marker = "!series_table_begin = Patient clinical parameters header descriptions"
    lines: list[str] = []
    active = False
    with gzip.open(path, "rt", encoding="utf-8", errors="strict") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line == marker:
                active = True
                continue
            if active and line == "!series_table_end":
                break
            if active:
                lines.append(line)
    if len(lines) != 287:
        raise ValueError(f"Expected header plus 286 clinical rows, found {len(lines)} lines")
    frame = pd.read_csv(io.StringIO("\n".join(lines)), sep="\t")
    frame = frame.rename(
        columns={
            "GEO asscession number": "gsm",
            "time to relapse or last follow-up (months)": "follow_up_months",
            "relapse (1=True)": "relapse",
            "ER Status": "er_status",
            "Brain relapses (1=yes, 0=no)": "brain_relapse",
            "lymph node status": "lymph_node_status",
        }
    )
    required = {"PID", "gsm", "follow_up_months", "relapse", "er_status"}
    if not required.issubset(frame.columns):
        raise ValueError("Clinical table columns changed")
    if frame["gsm"].duplicated().any() or len(frame) != 286:
        raise ValueError("Clinical sample identifiers are incomplete or duplicated")
    return frame


def construct_five_year_endpoint(clinical: pd.DataFrame) -> pd.DataFrame:
    """Create an auditable 60-month endpoint that respects early censoring.

    A relapse at or before 60 months is positive. A participant observed through
    60 months is negative at the 60-month landmark, even if relapse occurred
    later. Relapse-free participants censored before 60 months are excluded.
    """
    frame = clinical.copy()
    time = pd.to_numeric(frame["follow_up_months"], errors="raise")
    relapse = pd.to_numeric(frame["relapse"], errors="raise").astype(int)
    if not set(relapse.unique()).issubset({0, 1}) or (time <= 0).any():
        raise ValueError("Invalid relapse or follow-up values")
    positive = (relapse == 1) & (time <= 60)
    adequate_follow_up = time >= 60
    eligible = positive | adequate_follow_up
    frame["endpoint_eligible"] = eligible
    frame["relapse_by_60_months"] = np.where(positive, 1, np.where(adequate_follow_up, 0, np.nan))
    frame["endpoint_reason"] = np.select(
        [positive, adequate_follow_up, (relapse == 0) & (time <= 60)],
        ["relapse_by_60_months", "event_free_at_60_months", "censored_by_60_months"],
        default="invalid",
    )
    return frame


def parse_platform_annotation(path: Path) -> pd.DataFrame:
    lines: list[str] = []
    active = False
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line == "!platform_table_begin":
                active = True
                continue
            if active and line == "!platform_table_end":
                break
            if active:
                lines.append(line)
    if not lines:
        raise ValueError("GPL96 annotation table was not found")
    frame = pd.read_csv(io.StringIO("\n".join(lines)), sep="\t", dtype=str)
    required = {"ID", "Gene ID", "Gene symbol"}
    if not required.issubset(frame.columns):
        raise ValueError("GPL96 annotation columns changed")
    return frame


def collapse_to_unambiguous_genes(
    probes: pd.DataFrame, annotation: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Log-transform and select one outcome-blind representative probe per gene."""
    ann = annotation[["ID", "Gene ID", "Gene symbol"]].copy()
    ann.columns = ["probe_id", "entrez_id", "gene_symbol"]
    single = ann[ann["entrez_id"].str.fullmatch(r"\d+", na=False)].copy()
    single = single[~single["gene_symbol"].fillna("").str.contains("///", regex=False)]
    single = single[single["gene_symbol"].fillna("").str.strip().ne("")]
    single = single[single["probe_id"].isin(probes.index)].drop_duplicates("probe_id")
    values = np.log2(probes.loc[single["probe_id"]].astype(float) + 1.0)
    q75 = values.quantile(0.75, axis=1)
    q25 = values.quantile(0.25, axis=1)
    iqr = (q75 - q25).rename("iqr")
    selection = single.merge(iqr, left_on="probe_id", right_index=True)
    selection = selection.sort_values(["entrez_id", "iqr", "probe_id"], ascending=[True, False, True])
    selection = selection.drop_duplicates("entrez_id", keep="first")
    selection = selection.sort_values(["gene_symbol", "entrez_id", "probe_id"]).reset_index(drop=True)
    chosen = values.loc[selection["probe_id"]].copy()
    chosen.index = selection["gene_symbol"] + "|" + selection["entrez_id"]
    gene_by_sample = chosen
    if gene_by_sample.index.has_duplicates:
        # Symbols may map to multiple Entrez IDs, but the combined label is unique.
        raise ValueError("Combined gene labels are not unique")
    return gene_by_sample, selection


def validate_model_inputs(X: pd.DataFrame, y: pd.Series, sample_ids: Iterable[str]) -> None:
    ids = list(sample_ids)
    if len(X) != len(y) or len(X) != len(ids):
        raise ValueError("Expression, outcome, and sample identifiers differ in length")
    if len(set(ids)) != len(ids):
        raise ValueError("Sample identifiers are duplicated")
    if set(pd.Series(y).dropna().astype(int).unique()) != {0, 1}:
        raise ValueError("Outcome must contain both binary classes")
    if X.columns.duplicated().any() or X.index.duplicated().any():
        raise ValueError("Expression labels are duplicated")
    if not np.isfinite(X.to_numpy(dtype=float)).all():
        raise ValueError("Expression values must be finite")


def _expression_pipeline() -> Pipeline:
    return Pipeline(
        steps=[
            ("select", SelectKBest(score_func=f_classif)),
            ("scale", StandardScaler()),
            (
                "classifier",
                LogisticRegression(
                    solver="liblinear",
                    class_weight="balanced",
                    max_iter=3000,
                    random_state=RANDOM_SEED,
                ),
            ),
        ]
    )


def nested_repeated_cv(
    X: pd.DataFrame,
    y: pd.Series,
    sample_ids: list[str],
    outer_splits: int = 5,
    outer_repeats: int = 3,
    inner_splits: int = 4,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Run leakage-controlled nested CV and return predictions, folds, features."""
    validate_model_inputs(X, y, sample_ids)
    y_array = np.asarray(y, dtype=int)
    outer = RepeatedStratifiedKFold(
        n_splits=outer_splits, n_repeats=outer_repeats, random_state=RANDOM_SEED
    )
    param_grid = {
        "select__k": [25, 50, 100],
        "classifier__C": [0.01, 0.1, 1.0],
        "classifier__l1_ratio": [0.0, 1.0],
    }
    predictions: list[dict[str, object]] = []
    folds: list[dict[str, object]] = []
    features: list[dict[str, object]] = []
    n_splits_total = outer.get_n_splits(X, y_array)
    for fold_index, (train_idx, test_idx) in enumerate(outer.split(X, y_array), start=1):
        repeat = (fold_index - 1) // outer_splits + 1
        fold = (fold_index - 1) % outer_splits + 1
        inner = StratifiedKFold(n_splits=inner_splits, shuffle=True, random_state=RANDOM_SEED + fold_index)
        search = GridSearchCV(
            _expression_pipeline(),
            param_grid=param_grid,
            scoring="roc_auc",
            cv=inner,
            n_jobs=1,
            refit=True,
            error_score="raise",
        )
        search.fit(X.iloc[train_idx], y_array[train_idx])
        probability = search.predict_proba(X.iloc[test_idx])[:, 1]
        selected_mask = search.best_estimator_.named_steps["select"].get_support()
        selected_names = X.columns[selected_mask]
        coefficients = search.best_estimator_.named_steps["classifier"].coef_[0]
        for idx, prob in zip(test_idx, probability):
            predictions.append(
                {
                    "sample_id": sample_ids[idx],
                    "outcome": int(y_array[idx]),
                    "repeat": repeat,
                    "fold": fold,
                    "predicted_probability": float(prob),
                }
            )
        folds.append(
            {
                "repeat": repeat,
                "fold": fold,
                "train_n": len(train_idx),
                "test_n": len(test_idx),
                "test_events": int(y_array[test_idx].sum()),
                "inner_best_auc": float(search.best_score_),
                "test_auc": float(roc_auc_score(y_array[test_idx], probability)),
                "best_k": int(search.best_params_["select__k"]),
                "best_C": float(search.best_params_["classifier__C"]),
                "best_penalty": "l1" if search.best_params_["classifier__l1_ratio"] == 1.0 else "l2",
            }
        )
        for name, coef in zip(selected_names, coefficients):
            symbol, entrez = name.rsplit("|", 1)
            features.append(
                {
                    "repeat": repeat,
                    "fold": fold,
                    "feature": name,
                    "gene_symbol": symbol,
                    "entrez_id": entrez,
                    "coefficient": float(coef),
                }
            )
    pred = pd.DataFrame(predictions)
    if len(pred) != len(X) * outer_repeats or pred["sample_id"].value_counts().nunique() != 1:
        raise RuntimeError("Repeated-CV prediction coverage is incomplete")
    fold_frame = pd.DataFrame(folds)
    feature_frame = pd.DataFrame(features)
    if len(fold_frame) != n_splits_total:
        raise RuntimeError("Outer-fold accounting failed")
    return pred, fold_frame, feature_frame


def aggregate_repeated_predictions(predictions: pd.DataFrame) -> pd.DataFrame:
    grouped = (
        predictions.groupby(["sample_id", "outcome"], as_index=False)
        .agg(predicted_probability=("predicted_probability", "mean"),
             prediction_sd=("predicted_probability", "std"),
             n_predictions=("predicted_probability", "size"))
    )
    return grouped.sort_values("sample_id").reset_index(drop=True)


def _metrics(y: np.ndarray, probability: np.ndarray) -> dict[str, float]:
    predicted = (probability >= 0.5).astype(int)
    tn, fp, fn, tp = confusion_matrix(y, predicted, labels=[0, 1]).ravel()
    clipped = np.clip(probability, 1e-6, 1 - 1e-6)
    calibration = LogisticRegression(C=1e6, solver="lbfgs", max_iter=3000)
    calibration.fit(logit(clipped).reshape(-1, 1), y)
    return {
        "roc_auc": float(roc_auc_score(y, probability)),
        "average_precision": float(average_precision_score(y, probability)),
        "brier_score": float(brier_score_loss(y, probability)),
        "balanced_accuracy_0.5": float(balanced_accuracy_score(y, predicted)),
        "sensitivity_0.5": float(tp / (tp + fn)),
        "specificity_0.5": float(tn / (tn + fp)),
        "calibration_intercept": float(calibration.intercept_[0]),
        "calibration_slope": float(calibration.coef_[0, 0]),
    }


def bootstrap_metrics(
    y: np.ndarray, probability: np.ndarray, iterations: int = 2000
) -> pd.DataFrame:
    """Stratified participant bootstrap intervals for discrimination/calibration."""
    y = np.asarray(y, dtype=int)
    probability = np.asarray(probability, dtype=float)
    point = _metrics(y, probability)
    rng = np.random.default_rng(RANDOM_SEED)
    class_indices = [np.flatnonzero(y == value) for value in (0, 1)]
    draws: dict[str, list[float]] = {key: [] for key in point}
    for _ in range(iterations):
        sampled = np.concatenate([rng.choice(idx, size=len(idx), replace=True) for idx in class_indices])
        values = _metrics(y[sampled], probability[sampled])
        for key, value in values.items():
            draws[key].append(value)
    rows = []
    for metric, value in point.items():
        lower, upper = np.quantile(draws[metric], [0.025, 0.975])
        rows.append({"metric": metric, "estimate": value, "ci_lower": lower, "ci_upper": upper})
    return pd.DataFrame(rows)


def cross_validated_er_baseline(
    er_binary: np.ndarray, y: np.ndarray, sample_ids: list[str],
    outer_splits: int = 5, outer_repeats: int = 3,
) -> pd.DataFrame:
    er_binary = np.asarray(er_binary, dtype=float).reshape(-1, 1)
    y = np.asarray(y, dtype=int)
    outer = RepeatedStratifiedKFold(
        n_splits=outer_splits, n_repeats=outer_repeats, random_state=RANDOM_SEED
    )
    rows = []
    for fold_index, (train_idx, test_idx) in enumerate(outer.split(er_binary, y), start=1):
        repeat = (fold_index - 1) // outer_splits + 1
        fold = (fold_index - 1) % outer_splits + 1
        model = LogisticRegression(C=1e6, solver="liblinear", max_iter=3000)
        model.fit(er_binary[train_idx], y[train_idx])
        probs = model.predict_proba(er_binary[test_idx])[:, 1]
        for idx, prob in zip(test_idx, probs):
            rows.append({"sample_id": sample_ids[idx], "outcome": int(y[idx]),
                         "repeat": repeat, "fold": fold,
                         "predicted_probability": float(prob)})
    return pd.DataFrame(rows)


def cross_validated_prevalence_baseline(
    y: np.ndarray, sample_ids: list[str], outer_splits: int = 5, outer_repeats: int = 3
) -> pd.DataFrame:
    """Predict each held-out sample using only its training fold's event rate."""
    y = np.asarray(y, dtype=int)
    placeholder = np.zeros((len(y), 1), dtype=float)
    outer = RepeatedStratifiedKFold(
        n_splits=outer_splits, n_repeats=outer_repeats, random_state=RANDOM_SEED
    )
    rows = []
    for fold_index, (train_idx, test_idx) in enumerate(outer.split(placeholder, y), start=1):
        repeat = (fold_index - 1) // outer_splits + 1
        fold = (fold_index - 1) % outer_splits + 1
        probability = float(y[train_idx].mean())
        for idx in test_idx:
            rows.append({"sample_id": sample_ids[idx], "outcome": int(y[idx]),
                         "repeat": repeat, "fold": fold,
                         "predicted_probability": probability})
    return pd.DataFrame(rows)


def feature_stability(feature_rows: pd.DataFrame, total_folds: int) -> pd.DataFrame:
    if feature_rows.empty:
        return pd.DataFrame(columns=["feature", "gene_symbol", "entrez_id", "selection_count",
                                     "selection_fraction", "mean_coefficient", "sign_consistency"])
    grouped = feature_rows.groupby(["feature", "gene_symbol", "entrez_id"], as_index=False)
    result = grouped.agg(
        selection_count=("coefficient", "size"),
        mean_coefficient=("coefficient", "mean"),
        positive_fraction=("coefficient", lambda x: float((x > 0).mean())),
    )
    result["selection_fraction"] = result["selection_count"] / total_folds
    result["sign_consistency"] = np.maximum(result["positive_fraction"], 1 - result["positive_fraction"])
    return result.sort_values(["selection_fraction", "sign_consistency", "feature"],
                              ascending=[False, False, True]).reset_index(drop=True)


def curve_data(y: np.ndarray, probability: np.ndarray) -> tuple[pd.DataFrame, pd.DataFrame]:
    fpr, tpr, thresholds = roc_curve(y, probability)
    roc = pd.DataFrame({"false_positive_rate": fpr, "true_positive_rate": tpr,
                        "threshold": thresholds})
    observed, predicted = calibration_curve(y, probability, n_bins=8, strategy="quantile")
    calibration = pd.DataFrame({"mean_predicted_probability": predicted,
                                "observed_event_fraction": observed})
    return roc, calibration
