#!/usr/bin/env python3
"""Create a compact publication-oriented direct-interaction benchmark figure."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--output-prefix", required=True, type=Path)
    args = parser.parse_args()
    data = pd.read_csv(args.summary, sep="\t")
    data = data[(data["region_scope"] == "3UTR") & (data["rank_metric"] == "apa_exposure")].copy()
    order = ["hsa-miR-16-5p", "hsa-miR-17-5p", "hsa-miR-25-3p", "hsa-miR-92a-3p", "hsa-miR-186-5p"]
    labels = [x.replace("hsa-", "") for x in order]
    datasets = [
        ("GSE73057_CLEAR_CLIP", "CLEAR-CLIP (Huh-7.5)", "#0072B2", "o"),
        ("GSE50452_CLASH", "CLASH (HEK293)", "#D55E00", "s"),
    ]
    plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.5), constrained_layout=True)
    x = list(range(len(order)))
    for dataset, label, color, marker in datasets:
        subset = data[data.dataset == dataset].set_index("mirna").reindex(order)
        axes[0].plot(x, subset["positive_unlabeled_auc"], marker=marker, color=color, label=label, linewidth=1.5)
        axes[1].plot(x, subset["top_10_enrichment"], marker=marker, color=color, label=label, linewidth=1.5)
    axes[0].axhline(0.5, color="#777777", linestyle="--", linewidth=1)
    axes[1].axhline(1.0, color="#777777", linestyle="--", linewidth=1)
    axes[0].set_ylabel("Positive–unlabeled AUC")
    axes[1].set_ylabel("Top-10% enrichment (fold)")
    axes[0].set_ylim(0.4, 0.75)
    axes[1].set_ylim(0, 3.4)
    for ax, title in zip(axes, ["A  Rank distribution", "B  Top-decile recovery"]):
        ax.set_title(title, loc="left", fontweight="bold")
        ax.set_xticks(x, labels, rotation=35, ha="right")
        ax.grid(axis="y", color="#DDDDDD", linewidth=0.6)
    axes[0].legend(frameon=False, loc="lower left", fontsize=8)
    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output_prefix.with_suffix(".png"), dpi=300)
    fig.savefig(args.output_prefix.with_suffix(".svg"))
    fig.savefig(args.output_prefix.with_suffix(".pdf"))


if __name__ == "__main__":
    main()
