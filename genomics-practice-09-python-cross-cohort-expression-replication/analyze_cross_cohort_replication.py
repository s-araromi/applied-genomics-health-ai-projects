#!/usr/bin/env python3
"""Cross-platform replication of paired lung-tumor expression differences."""

from __future__ import annotations

import argparse
import json
import platform
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
from scipy import stats

from expression_replication import (
    SOURCES,
    benjamini_hochberg,
    clustered_ols,
    complete_pair_indices,
    ensure_sources,
    fixed_effect_meta,
    paired_statistics,
    read_annotation,
    read_geo_matrix,
    representative_genes,
    validate_output_accounting,
)

DISCOVERY_Q = 0.05
DISCOVERY_EFFECT = 1.0
VALIDATION_Q = 0.05
VALIDATION_EFFECT = 0.5


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, default=Path("data/cache"))
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"))
    parser.add_argument("--no-download", action="store_true",
                        help="Fail if a pinned source is absent instead of downloading it")
    return parser.parse_args()


def build_results(cache: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, object]]:
    discovery_x, discovery_m = read_geo_matrix(
        cache / "GSE19804_series_matrix.txt.gz", series="GSE19804", platform="GPL570",
        expected_shape=(54_675, 120),
    )
    validation_x, validation_m = read_geo_matrix(
        cache / "GSE10072_series_matrix.txt.gz", series="GSE10072", platform="GPL96",
        expected_shape=(22_283, 107),
    )
    discovery_a = read_annotation(cache / "GPL570.annot.gz", 54_675)
    validation_a = read_annotation(cache / "GPL96.annot.gz", 22_283)
    discovery_r = representative_genes(discovery_x, discovery_a).set_index("gene_id")
    validation_r = representative_genes(validation_x, validation_a).set_index("gene_id")
    common = sorted(set(discovery_r.index) & set(validation_r.index), key=int)
    if len(common) < 8_000:
        raise ValueError("Unexpectedly small cross-platform gene intersection")

    d_probe = discovery_r.loc[common, "probe_id"]
    v_probe = validation_r.loc[common, "probe_id"]
    d_values = discovery_x.loc[d_probe].to_numpy(dtype=float)
    v_values = validation_x.loc[v_probe].to_numpy(dtype=float)
    d_tumor, d_normal, d_patients = complete_pair_indices(discovery_m)
    v_tumor, v_normal, v_patients = complete_pair_indices(validation_m)
    if len(d_patients) != 60 or len(v_patients) != 33:
        raise ValueError("Unexpected complete-pair counts")

    d = paired_statistics(d_values, d_tumor, d_normal).add_prefix("discovery_")
    v = paired_statistics(v_values, v_tumor, v_normal).add_prefix("validation_")
    sensitivity = clustered_ols(v_values, validation_m)
    results = pd.DataFrame({
        "gene_id": common,
        "gene_symbol": discovery_r.loc[common, "gene_symbol"].to_numpy(),
        "gene_title": discovery_r.loc[common, "gene_title"].to_numpy(),
        "discovery_probe": d_probe.to_numpy(),
        "validation_probe": v_probe.to_numpy(),
        "discovery_probe_iqr": discovery_r.loc[common, "probe_iqr"].to_numpy(),
        "validation_probe_iqr": validation_r.loc[common, "probe_iqr"].to_numpy(),
    })
    results = pd.concat([results, d, v, sensitivity], axis=1)
    results["discovery_q"] = benjamini_hochberg(results["discovery_p_value"].to_numpy())
    results["discovery_candidate"] = (
        (results["discovery_q"] <= DISCOVERY_Q)
        & (results["discovery_mean_log2_difference"].abs() >= DISCOVERY_EFFECT)
    )
    results["validation_q_within_candidates"] = 1.0
    candidates = results["discovery_candidate"].to_numpy()
    results.loc[candidates, "validation_q_within_candidates"] = benjamini_hochberg(
        results.loc[candidates, "validation_p_value"].to_numpy()
    )
    same_direction = np.sign(results["discovery_mean_log2_difference"]) == np.sign(
        results["validation_mean_log2_difference"]
    )
    results["validation_same_direction"] = same_direction
    results["clustered_sensitivity_pass"] = (
        (np.sign(results["discovery_mean_log2_difference"]) == np.sign(results["clustered_effect"]))
        & (results["clustered_p_value"] < 0.05)
    )
    results["replicated"] = (
        results["discovery_candidate"]
        & results["validation_same_direction"]
        & (results["validation_q_within_candidates"] <= VALIDATION_Q)
        & (results["validation_mean_log2_difference"].abs() >= VALIDATION_EFFECT)
        & results["clustered_sensitivity_pass"]
    )
    validate_output_accounting(results)

    replicated = results.loc[results["replicated"]].copy()
    meta = fixed_effect_meta(
        replicated["discovery_mean_log2_difference"].to_numpy(),
        replicated["discovery_standard_error"].to_numpy(),
        replicated["validation_mean_log2_difference"].to_numpy(),
        replicated["validation_standard_error"].to_numpy(),
    )
    replicated = pd.concat([replicated.reset_index(drop=True), meta], axis=1)
    replicated["absolute_fixed_effect"] = replicated["fixed_effect_log2_difference"].abs()
    replicated = replicated.sort_values(
        ["absolute_fixed_effect", "fixed_effect_p_value", "gene_id"],
        ascending=[False, True, True],
    ).drop(columns="absolute_fixed_effect")

    cohort_samples = pd.concat([
        discovery_m.assign(cohort="GSE19804", primary_pair=discovery_m.patient_id.isin(d_patients)),
        validation_m.assign(cohort="GSE10072", primary_pair=validation_m.patient_id.isin(v_patients)),
    ], ignore_index=True)
    candidate_effects = results.loc[candidates, [
        "discovery_mean_log2_difference", "validation_mean_log2_difference"
    ]]
    rho, rho_p = stats.spearmanr(
        results["discovery_mean_log2_difference"], results["validation_mean_log2_difference"]
    )
    summary = {
        "common_entrez_genes": len(results),
        "discovery_complete_pairs": len(d_patients),
        "validation_complete_pairs": len(v_patients),
        "validation_all_samples": len(validation_m),
        "validation_participants": validation_m.patient_id.nunique(),
        "discovery_candidates": int(candidates.sum()),
        "candidate_same_direction": int((same_direction & candidates).sum()),
        "replicated_genes": len(replicated),
        "replication_fraction": len(replicated) / int(candidates.sum()),
        "all_gene_spearman_rho": float(rho),
        "all_gene_spearman_p": float(rho_p),
        "candidate_pearson_r": float(candidate_effects.corr(method="pearson").iloc[0, 1]),
    }
    return results, replicated, cohort_samples, summary


def plot_results(results: pd.DataFrame, replicated: pd.DataFrame, output_dir: Path) -> None:
    candidates = results.loc[results["discovery_candidate"]].copy()
    figure, axes = plt.subplots(1, 2, figsize=(13.5, 6.4), constrained_layout=True)
    ax = axes[0]
    other = candidates.loc[~candidates["replicated"]]
    hit = candidates.loc[candidates["replicated"]]
    ax.scatter(other["discovery_mean_log2_difference"], other["validation_mean_log2_difference"],
               s=18, alpha=0.45, color="#778899", edgecolor="none", label="Not replicated")
    ax.scatter(hit["discovery_mean_log2_difference"], hit["validation_mean_log2_difference"],
               s=25, alpha=0.8, color="#B4493D", edgecolor="white", linewidth=0.3,
               label="Replicated")
    limits = np.array([
        candidates[["discovery_mean_log2_difference", "validation_mean_log2_difference"]].min().min(),
        candidates[["discovery_mean_log2_difference", "validation_mean_log2_difference"]].max().max(),
    ])
    pad = 0.08 * (limits[1] - limits[0])
    limits += np.array([-pad, pad])
    ax.plot(limits, limits, linestyle="--", color="#333333", linewidth=1, label="Equal effect")
    ax.axhline(0, color="#BBBBBB", linewidth=0.8)
    ax.axvline(0, color="#BBBBBB", linewidth=0.8)
    ax.set_xlim(limits); ax.set_ylim(limits)
    ax.set_xlabel("GSE19804 paired mean log2 difference")
    ax.set_ylabel("GSE10072 paired mean log2 difference")
    ax.set_title("Discovery candidates across cohorts")
    ax.legend(frameon=False, loc="upper left")

    ax = axes[1]
    top = replicated.head(12).sort_values("fixed_effect_log2_difference")
    y = np.arange(len(top))
    ax.errorbar(top["fixed_effect_log2_difference"], y,
                xerr=1.96 * top["fixed_effect_standard_error"], fmt="o",
                color="#243B53", ecolor="#829AB1", capsize=3, label="Fixed effect")
    ax.scatter(top["discovery_mean_log2_difference"], y - 0.15, marker="|", s=80,
               color="#2A7F62", label="GSE19804")
    ax.scatter(top["validation_mean_log2_difference"], y + 0.15, marker="|", s=80,
               color="#D98732", label="GSE10072")
    ax.axvline(0, color="#777777", linewidth=0.8)
    ax.set_yticks(y, top["gene_symbol"])
    ax.set_xlabel("Tumor − normal log2 expression")
    ax.set_title("Largest replicated fixed effects")
    ax.legend(frameon=False, loc="best")
    figure.suptitle("Cross-platform replication of lung-tumor expression differences", fontsize=14)
    for suffix in ("png", "svg"):
        path = output_dir / f"cross_cohort_replication.{suffix}"
        figure.savefig(path, dpi=180 if suffix == "png" else None)
        if suffix == "svg":
            # Matplotlib's SVG serializer adds trailing spaces to path lines.
            # Normalizing them keeps Git's whitespace validation meaningful.
            lines = path.read_text(encoding="utf-8").splitlines()
            path.write_text("\n".join(line.rstrip() for line in lines) + "\n", encoding="utf-8")
    plt.close(figure)


def write_outputs(output_dir: Path, results: pd.DataFrame, replicated: pd.DataFrame,
                  cohort_samples: pd.DataFrame, summary: dict[str, object]) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    results.to_csv(output_dir / "cross_cohort_gene_results.csv", index=False, lineterminator="\n")
    replicated.to_csv(output_dir / "replicated_gene_results.csv", index=False, lineterminator="\n")
    cohort_samples.to_csv(output_dir / "cohort_sample_manifest.csv", index=False, lineterminator="\n")
    summary_values = [str(value) if isinstance(value, int) else format(value, ".15g")
                      for value in summary.values()]
    pd.DataFrame({"metric": summary.keys(), "value": summary_values}).to_csv(
        output_dir / "analysis_summary.csv", index=False, lineterminator="\n"
    )
    top = replicated.head(10)
    lines = [
        "# Project 09: Cross-Cohort Expression Replication", "",
        "## Verified analysis accounting", "",
        f"- Common Entrez Gene IDs tested: {summary['common_entrez_genes']:,}",
        f"- Discovery pairs (GSE19804): {summary['discovery_complete_pairs']}",
        f"- Validation pairs (GSE10072): {summary['validation_complete_pairs']}",
        f"- Discovery candidates: {summary['discovery_candidates']:,}",
        f"- Replicated genes: {summary['replicated_genes']:,} ({100 * summary['replication_fraction']:.1f}%)", "",
        "## Cross-cohort agreement", "",
        f"All-gene Spearman rho: {summary['all_gene_spearman_rho']:.4f}.",
        f"Discovery-candidate Pearson r: {summary['candidate_pearson_r']:.4f}.", "",
        "## Largest replicated fixed effects", "",
        "| Gene | GSE19804 effect | GSE10072 effect | Fixed effect | I2 |",
        "| --- | ---: | ---: | ---: | ---: |",
    ]
    for row in top.itertuples():
        lines.append(
            f"| {row.gene_symbol} | {row.discovery_mean_log2_difference:.3f} | "
            f"{row.validation_mean_log2_difference:.3f} | {row.fixed_effect_log2_difference:.3f} | "
            f"{row.heterogeneity_I2:.3f} |"
        )
    lines.extend([
        "", "## Interpretation boundary", "",
        "Replication means that a pre-specified statistical and effect-size rule was met in this second cohort.",
        "It does not establish diagnostic performance, causality, clinical utility, or independence from tissue composition.",
        "The fixed-effect summaries combine two platforms and populations and should be read with the reported heterogeneity.",
    ])
    (output_dir / "replication_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    metadata = {
        "python": sys.version.split()[0], "platform": platform.platform(),
        "numpy": np.__version__, "pandas": pd.__version__,
        "scipy": scipy.__version__, "matplotlib": matplotlib.__version__,
        "thresholds": {
            "discovery_q": DISCOVERY_Q, "discovery_abs_effect": DISCOVERY_EFFECT,
            "validation_q_within_candidates": VALIDATION_Q,
            "validation_abs_effect": VALIDATION_EFFECT,
            "clustered_sensitivity_p": 0.05,
        },
    }
    (output_dir / "run_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")


def main() -> None:
    args = parse_args()
    ensure_sources(args.cache, allow_download=not args.no_download)
    results, replicated, cohort_samples, summary = build_results(args.cache)
    write_outputs(args.output_dir, results, replicated, cohort_samples, summary)
    plot_results(results, replicated, args.output_dir)
    print(pd.Series(summary).to_string())


if __name__ == "__main__":
    main()
