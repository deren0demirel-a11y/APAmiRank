#!/usr/bin/env python3
"""Create the publication figure for the locked luciferase holdout."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


APA = "consensus_apa_exposure_percentile"
UNWEIGHTED = "consensus_unweighted_additive_percentile"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pair-ranks", required=True, type=Path)
    parser.add_argument("--screening", required=True, type=Path)
    parser.add_argument("--null-draws", required=True, type=Path)
    parser.add_argument("--uncertainty-summary", required=True, type=Path)
    parser.add_argument("--output-prefix", required=True, type=Path)
    args = parser.parse_args()

    pairs = pd.read_csv(args.pair_ranks, sep="\t")
    screening = pd.read_csv(args.screening, sep="\t")
    observed = pairs.dropna(subset=[APA, UNWEIGHTED]).copy().sort_values(APA)
    null = pd.read_csv(args.null_draws, sep="\t")
    uncertainty = pd.read_csv(args.uncertainty_summary, sep="\t").set_index("metric")
    blue, orange, grey = "#0072B2", "#D55E00", "#666666"
    plt.rcParams.update({
        "font.size": 8.5,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "pdf.fonttype": 42,
        "svg.fonttype": "none",
        "svg.hashsalt": "APAmiRank-v0.2.18",
    })
    fig, axes = plt.subplots(2, 2, figsize=(8.4, 7.0), constrained_layout=True)

    ax = axes[0, 0]
    y = np.arange(len(observed))
    ax.scatter(observed[APA], y, color=blue, s=28, zorder=3)
    ax.axvline(0.5, color=grey, linestyle="--", linewidth=0.9)
    ax.axvline(0.8, color="#999999", linestyle=":", linewidth=0.8)
    ax.axvline(0.9, color="#999999", linestyle=":", linewidth=0.8)
    ax.set_yticks(y, observed.gene)
    ax.set_xlim(0, 1.02)
    ax.set_xlabel("APA-exposure percentile")
    ax.set_title("A  Independently screened targets", loc="left", fontweight="bold")
    ax.grid(axis="x", color="#E2E2E2", linewidth=0.6)

    ax = axes[0, 1]
    ax.hist(null.mean_apa, bins=45, color="#B9D8E8", edgecolor="white", linewidth=0.3)
    observed_mean = observed[APA].mean()
    ax.axvline(observed_mean, color=orange, linewidth=2, label=f"Observed = {observed_mean:.3f}")
    ax.axvline(null.mean_apa.mean(), color=grey, linestyle="--", linewidth=1, label="Null mean")
    adjusted_p = uncertainty.loc["mean_apa", "matched_null_holm_p"]
    ax.text(0.03, 0.94, f"Holm-adjusted p = {adjusted_p:.4f}", transform=ax.transAxes,
            ha="left", va="top", fontsize=7.5)
    ax.set_xlabel("Mean APA-exposure percentile")
    ax.set_ylabel("Monte Carlo draws")
    ax.set_title("B  miRNA-matched rank-uniform null", loc="left", fontweight="bold")
    ax.legend(frameon=False, fontsize=7.5)

    ax = axes[1, 0]
    ax.scatter(observed[UNWEIGHTED], observed[APA], color=blue, s=30)
    ax.plot([0, 1], [0, 1], color=grey, linestyle="--", linewidth=0.9)
    offsets = {
        "CPEB2": (4, 8), "WEE1": (4, -2), "MDM2": (4, -11),
        "STAT3": (4, 3), "GRB2": (-10, 5), "PPM1D": (4, -8),
    }
    for row in observed.itertuples():
        offset = offsets.get(row.gene, (4, 3))
        ax.annotate(row.gene, (getattr(row, UNWEIGHTED), getattr(row, APA)),
                    xytext=offset, textcoords="offset points", fontsize=6.5)
    ax.set_xlim(0, 1.02)
    ax.set_ylim(0, 1.02)
    ax.set_xlabel("Unweighted-additive percentile")
    ax.set_ylabel("APA-exposure percentile")
    ax.set_title("C  Paired rank comparison", loc="left", fontweight="bold")
    ax.grid(color="#E2E2E2", linewidth=0.6)

    ax = axes[1, 1]
    labels = ["Locked", "Strictly included", "Observable"]
    counts = [len(screening), int((screening.eligibility_decision == "include").sum()), len(observed)]
    ax.bar(labels, counts, color=["#999999", blue, orange], width=0.65)
    for index, count in enumerate(counts):
        ax.text(index, count + 0.7, str(count), ha="center", fontweight="bold")
    ax.set_ylim(0, max(counts) * 1.16)
    ax.set_ylabel("Pairs")
    ax.set_title("D  Screening and rank coverage", loc="left", fontweight="bold")
    ax.tick_params(axis="x", rotation=15)
    ax.grid(axis="y", color="#E2E2E2", linewidth=0.6)

    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
    fixed_time = datetime(2026, 9, 10, tzinfo=timezone.utc)
    fig.savefig(
        args.output_prefix.with_suffix(".png"), dpi=300,
        metadata={"Software": "APAmiRank"},
    )
    fig.savefig(
        args.output_prefix.with_suffix(".svg"),
        metadata={"Creator": "APAmiRank", "Date": "2026-09-10"},
    )
    fig.savefig(
        args.output_prefix.with_suffix(".pdf"),
        metadata={"Creator": "APAmiRank", "CreationDate": fixed_time, "ModDate": fixed_time},
    )


if __name__ == "__main__":
    main()
