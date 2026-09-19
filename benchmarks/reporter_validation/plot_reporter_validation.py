#!/usr/bin/env python3
"""Plot the reporter benchmark and APA increment with bootstrap intervals."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ORDER = [
    "hsa-miR-16-5p", "hsa-miR-17-5p", "hsa-miR-25-3p",
    "hsa-miR-92a-3p", "hsa-miR-124-3p", "hsa-miR-155-5p",
    "hsa-miR-186-5p",
]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--increment", required=True, type=Path)
    parser.add_argument("--output-prefix", required=True, type=Path)
    args = parser.parse_args()
    summary = pd.read_csv(args.summary, sep="\t")
    primary = summary[
        (summary.universe_scope == "at_least_2_of_4_contexts")
        & (summary.evidence_set == "primary_reporter_plus_wb_or_qpcr")
        & (summary.rank_metric == "apa_exposure")
    ].set_index("miRNA").reindex(ORDER)
    increment = pd.read_csv(args.increment, sep="\t")
    increment = increment[
        (increment.universe_scope == "at_least_2_of_4_contexts")
        & (increment.evidence_set == "primary_reporter_plus_wb_or_qpcr")
    ].set_index("miRNA").reindex(ORDER)
    labels = [name.replace("hsa-", "") for name in ORDER]
    x = np.arange(len(ORDER))
    blue, orange = "#0072B2", "#D55E00"
    plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(1, 3, figsize=(11.0, 3.6), constrained_layout=True)
    axes[0].bar(x, primary.positive_unlabeled_auc, color=blue, width=0.68)
    axes[0].axhline(0.5, color="#666666", linestyle="--", linewidth=1)
    axes[0].set_ylim(0.5, 0.95)
    axes[0].set_ylabel("Positive–unlabeled AUC")
    for i, n in enumerate(primary.overlapping_positive_genes.astype(int)):
        axes[0].text(i, primary.positive_unlabeled_auc.iloc[i] + 0.012, f"n={n}", ha="center", fontsize=7)
    axes[1].bar(x, primary.top_10_enrichment, color=blue, width=0.68)
    axes[1].axhline(1.0, color="#666666", linestyle="--", linewidth=1)
    axes[1].set_ylim(0, 7.4)
    axes[1].set_ylabel("Top-10% enrichment (fold)")
    delta = increment.apa_minus_unweighted_auc.to_numpy(dtype=float)
    low = increment.bootstrap_ci_low.to_numpy(dtype=float)
    high = increment.bootstrap_ci_high.to_numpy(dtype=float)
    valid = np.isfinite(delta)
    axes[2].errorbar(x[valid], delta[valid], yerr=[delta[valid] - low[valid], high[valid] - delta[valid]], fmt="o", color=orange, ecolor=orange, capsize=3)
    axes[2].axhline(0, color="#666666", linestyle="--", linewidth=1)
    axes[2].set_ylim(-0.09, 0.15)
    axes[2].set_ylabel("APA minus unweighted PU-AUC")
    axes[2].text(x[-1], 0.01, "n<5", ha="center", va="bottom", fontsize=7, color="#555555")
    for ax, title in zip(axes, ["A  Reporter recovery", "B  Top-decile enrichment", "C  Increment from APA exposure"]):
        ax.set_title(title, loc="left", fontweight="bold")
        ax.set_xticks(x, labels, rotation=40, ha="right")
        ax.grid(axis="y", color="#DDDDDD", linewidth=0.6)
    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output_prefix.with_suffix(".png"), dpi=300)
    fig.savefig(args.output_prefix.with_suffix(".svg"))
    fig.savefig(args.output_prefix.with_suffix(".pdf"))


if __name__ == "__main__":
    main()
