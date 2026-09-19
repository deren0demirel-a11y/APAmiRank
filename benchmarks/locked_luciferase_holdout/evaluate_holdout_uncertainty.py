#!/usr/bin/env python3
"""Rank-uniform matched-null tests and miRNA-cluster bootstrap for the holdout."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


APA = "consensus_apa_exposure_percentile"
UNWEIGHTED = "consensus_unweighted_additive_percentile"
METRICS = ["mean_apa", "median_apa", "top_decile_fraction", "top_quintile_fraction"]


def holdout_statistics(apa: np.ndarray) -> dict[str, float]:
    return {
        "mean_apa": float(np.mean(apa)),
        "median_apa": float(np.median(apa)),
        "top_decile_fraction": float(np.mean(apa >= 0.90)),
        "top_quintile_fraction": float(np.mean(apa >= 0.80)),
    }


def holm_adjust(p_values: pd.Series) -> pd.Series:
    valid = p_values.dropna().sort_values()
    adjusted = pd.Series(np.nan, index=p_values.index, dtype=float)
    if valid.empty:
        return adjusted
    scaled = valid.to_numpy() * np.arange(len(valid), 0, -1)
    scaled = np.maximum.accumulate(scaled).clip(max=1.0)
    adjusted.loc[valid.index] = scaled
    return adjusted


def simulate_rank_uniform_null(
    cluster_counts: dict[str, int],
    universe_sizes: dict[str, int],
    iterations: int,
    seed: int,
) -> pd.DataFrame:
    """Sample rank positions without replacement within each miRNA universe."""
    rng = np.random.default_rng(seed)
    rows = []
    for iteration in range(iterations):
        sampled = []
        for mirna, count in cluster_counts.items():
            size = universe_sizes[mirna]
            if count > size or size < 2:
                raise ValueError(f"invalid matched-null size for {mirna}: k={count}, N={size}")
            zero_based_ranks = rng.choice(size, size=count, replace=False)
            sampled.append(1.0 - zero_based_ranks / (size - 1.0))
        row = holdout_statistics(np.concatenate(sampled))
        row["iteration"] = iteration + 1
        rows.append(row)
    return pd.DataFrame(rows)[["iteration", *METRICS]]


def cluster_table(observed: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for mirna, group in observed.groupby("miRNA", sort=True):
        apa = group[APA].to_numpy(dtype=float)
        delta = (group[APA] - group[UNWEIGHTED]).to_numpy(dtype=float)
        row = {"miRNA": mirna, "observable_pairs": len(group), **holdout_statistics(apa)}
        row["median_apa_minus_unweighted"] = float(np.median(delta))
        rows.append(row)
    return pd.DataFrame(rows)


def bootstrap_clusters(clusters: pd.DataFrame, iterations: int, seed: int) -> pd.DataFrame:
    """Resample complete miRNA clusters and give every sampled cluster equal weight."""
    metric_columns = [*METRICS, "median_apa_minus_unweighted"]
    values = clusters[metric_columns].to_numpy(dtype=float)
    rng = np.random.default_rng(seed)
    indices = rng.integers(0, len(values), size=(iterations, len(values)))
    draws = values[indices].mean(axis=1)
    return pd.DataFrame(draws, columns=metric_columns).assign(iteration=np.arange(1, iterations + 1))[
        ["iteration", *metric_columns]
    ]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pair-ranks", required=True, type=Path)
    parser.add_argument("--universe-summary", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--iterations", type=int, default=50_000)
    parser.add_argument("--seed", type=int, default=20_260_910)
    args = parser.parse_args()

    pairs = pd.read_csv(args.pair_ranks, sep="\t")
    observed = pairs.dropna(subset=[APA, UNWEIGHTED]).copy()
    if observed.empty:
        raise ValueError("no observable pairs with both APA and unweighted percentiles")
    if observed[["miRNA", "gene"]].duplicated().any():
        raise ValueError("observable holdout pairs must be unique")

    universe = pd.read_csv(args.universe_summary, sep="\t")
    universe = universe[
        universe.universe_scope.eq("at_least_2_of_4_contexts")
        & universe.evidence_set.eq("primary_reporter_plus_wb_or_qpcr")
        & universe.rank_metric.eq("apa_exposure")
    ]
    if universe.miRNA.duplicated().any():
        raise ValueError("universe summary contains duplicate primary miRNA rows")
    universe_sizes = universe.set_index("miRNA").candidate_genes.astype(int).to_dict()
    cluster_counts = observed.groupby("miRNA").size().astype(int).to_dict()
    missing = sorted(set(cluster_counts) - set(universe_sizes))
    if missing:
        raise ValueError(f"missing universe sizes for: {', '.join(missing)}")

    null = simulate_rank_uniform_null(cluster_counts, universe_sizes, args.iterations, args.seed)
    clusters = cluster_table(observed)
    bootstrap = bootstrap_clusters(clusters, args.iterations, args.seed + 1)
    observed_stats = holdout_statistics(observed[APA].to_numpy(dtype=float))

    rows = []
    for metric in METRICS:
        value = observed_stats[metric]
        p_value = (1 + int((null[metric] >= value).sum())) / (args.iterations + 1)
        rows.append({
            "metric": metric,
            "pooled_observed": value,
            "matched_null_mean": float(null[metric].mean()),
            "matched_null_one_sided_p": p_value,
            "cluster_equal_estimate": float(clusters[metric].mean()),
            "cluster_bootstrap_ci_low": float(bootstrap[metric].quantile(0.025)),
            "cluster_bootstrap_ci_high": float(bootstrap[metric].quantile(0.975)),
        })
    delta_metric = "median_apa_minus_unweighted"
    pooled_delta = float(np.median(observed[APA] - observed[UNWEIGHTED]))
    rows.append({
        "metric": delta_metric,
        "pooled_observed": pooled_delta,
        "matched_null_mean": np.nan,
        "matched_null_one_sided_p": np.nan,
        "cluster_equal_estimate": float(clusters[delta_metric].mean()),
        "cluster_bootstrap_ci_low": float(bootstrap[delta_metric].quantile(0.025)),
        "cluster_bootstrap_ci_high": float(bootstrap[delta_metric].quantile(0.975)),
    })
    summary = pd.DataFrame(rows)
    summary["matched_null_holm_p"] = holm_adjust(summary.matched_null_one_sided_p)
    summary["null_iterations"] = args.iterations
    summary["bootstrap_iterations"] = args.iterations
    summary["random_seed"] = args.seed
    summary["observable_pairs"] = len(observed)
    summary["miRNA_clusters"] = observed.miRNA.nunique()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    summary.to_csv(args.output_dir / "holdout_uncertainty_summary.tsv", sep="\t", index=False)
    clusters.to_csv(args.output_dir / "holdout_mirna_cluster_summary.tsv", sep="\t", index=False)
    deterministic_gzip = {"method": "gzip", "mtime": 0}
    null.to_csv(
        args.output_dir / "holdout_rank_uniform_null.tsv.gz",
        sep="\t", index=False, compression=deterministic_gzip,
    )
    bootstrap.to_csv(
        args.output_dir / "holdout_cluster_bootstrap.tsv.gz",
        sep="\t", index=False, compression=deterministic_gzip,
    )
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
