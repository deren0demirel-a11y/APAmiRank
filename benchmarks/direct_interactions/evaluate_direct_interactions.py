#!/usr/bin/env python3
"""Evaluate APAmiRank ranks against positive-only CLASH/CLEAR-CLIP evidence."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import hypergeom, mannwhitneyu


def percentile(rank: pd.Series, n: int) -> pd.Series:
    if n <= 1:
        return pd.Series(np.ones(len(rank)), index=rank.index)
    return 1.0 - (rank.astype(float) - 1.0) / (n - 1.0)


def evaluate_one(runs: Path, interactions: pd.DataFrame, dataset: str, cell: str, mirna: str, region_filter: str):
    run = runs / f"{cell}__{mirna}"
    base = pd.read_csv(run / "03_gene_ranking.tsv", sep="\t")
    cond = pd.read_csv(run / "03c_condition_gene_ranking.tsv", sep="\t")
    sites = pd.read_csv(run / "02_site_ranking.tsv", sep="\t")
    n = len(base)
    base["best_site_percentile"] = percentile(base["best_equal_weight_rank"], n)
    base["opportunity_percentile"] = percentile(base["opportunity_adjusted_descriptive_rank"], n)
    additive = sites.groupby("gene", as_index=False)["binding_evidence_score_equal_weight"].sum().rename(columns={"binding_evidence_score_equal_weight": "unweighted_additive_score"})
    additive["unweighted_additive_rank"] = additive["unweighted_additive_score"].rank(method="min", ascending=False).astype(int)
    additive["unweighted_additive_percentile"] = percentile(additive["unweighted_additive_rank"], n)
    cond = cond.rename(columns={"condition_gene_rank": "apa_exposure_rank"})
    cond["apa_exposure_percentile"] = percentile(cond["apa_exposure_rank"], n)
    ranks = base.merge(additive, on="gene", how="left").merge(cond[["gene", "apa_exposure_rank", "apa_exposure_percentile"]], on="gene", how="left")
    ranks["apa_percentile_gain"] = ranks["apa_exposure_percentile"] - ranks["unweighted_additive_percentile"]

    pos = interactions[(interactions.dataset == dataset) & (interactions.mirna == mirna)]
    if region_filter == "3UTR":
        pos = pos[pos.region.astype(str).str.contains("3'UTR", regex=False)]
    positive_genes = set(pos.gene)
    ranks["direct_positive"] = ranks.gene.isin(positive_genes)
    k = int(ranks.direct_positive.sum())
    m = n - k
    metrics = []
    for metric, col in [
        ("best_site", "best_site_percentile"),
        ("opportunity_adjusted", "opportunity_percentile"),
        ("unweighted_additive", "unweighted_additive_percentile"),
        ("apa_exposure", "apa_exposure_percentile"),
    ]:
        positives = ranks.loc[ranks.direct_positive, col].dropna()
        negatives = ranks.loc[~ranks.direct_positive, col].dropna()
        if len(positives) and len(negatives):
            u = mannwhitneyu(positives, negatives, alternative="greater")
            pu_auc = u.statistic / (len(positives) * len(negatives))
            p_rank = u.pvalue
        else:
            pu_auc = p_rank = np.nan
        row = {
            "dataset": dataset, "cell_line_rank": cell, "mirna": mirna,
            "region_scope": region_filter, "rank_metric": metric,
            "candidate_genes": n, "source_positive_genes": len(positive_genes),
            "overlapping_positive_genes": len(positives),
            "candidate_coverage_of_source_positives": len(positives) / len(positive_genes) if positive_genes else np.nan,
            "median_positive_percentile": positives.median() if len(positives) else np.nan,
            "positive_unlabeled_auc": pu_auc, "mann_whitney_one_sided_p": p_rank,
        }
        for frac in (0.10, 0.20):
            top_n = max(1, int(np.ceil(frac * n)))
            hits = int((positives >= (1 - (top_n - 1) / (n - 1))).sum()) if n > 1 else int(len(positives))
            row[f"top_{int(frac*100)}_hits"] = hits
            row[f"top_{int(frac*100)}_source_recall"] = hits / len(positive_genes) if positive_genes else np.nan
            row[f"top_{int(frac*100)}_observable_recall"] = hits / len(positives) if len(positives) else np.nan
            row[f"top_{int(frac*100)}_enrichment"] = (hits / len(positives)) / (top_n / n) if len(positives) else np.nan
            row[f"top_{int(frac*100)}_hypergeom_p"] = hypergeom.sf(hits - 1, n, len(positives), top_n) if len(positives) else np.nan
        metrics.append(row)
    detail = ranks[ranks.direct_positive].copy()
    detail.insert(0, "dataset", dataset)
    detail.insert(1, "region_scope", region_filter)
    detail.insert(2, "cell_line_rank", cell)
    detail.insert(3, "benchmark_mirna", mirna)
    return metrics, detail


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--interactions", required=True, type=Path)
    parser.add_argument("--runs", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    interactions = pd.read_csv(args.interactions, sep="\t")
    panel = ["hsa-miR-16-5p", "hsa-miR-17-5p", "hsa-miR-25-3p", "hsa-miR-92a-3p", "hsa-miR-186-5p"]
    specs = [("GSE73057_CLEAR_CLIP", "Huh7", m, "3UTR") for m in panel]
    specs += [("GSE50452_CLASH", "HEK293", m, "3UTR") for m in panel]
    specs += [("GSE50452_CLASH", "HEK293", m, "all_mRNA") for m in panel]
    rows, details = [], []
    for spec in specs:
        one_rows, one_detail = evaluate_one(args.runs, interactions, *spec)
        rows.extend(one_rows); details.append(one_detail)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    summary = pd.DataFrame(rows)
    summary["mann_whitney_bh_q"] = np.nan
    for _, idx in summary.groupby(["dataset", "region_scope", "rank_metric"]).groups.items():
        valid = summary.loc[idx, "mann_whitney_one_sided_p"].dropna().sort_values()
        if valid.empty:
            continue
        raw = valid.to_numpy() * len(valid) / np.arange(1, len(valid) + 1)
        adjusted = np.minimum.accumulate(raw[::-1])[::-1].clip(max=1.0)
        summary.loc[valid.index, "mann_whitney_bh_q"] = adjusted
    summary.to_csv(args.output_dir / "direct_interaction_benchmark_summary.tsv", sep="\t", index=False)
    pd.concat(details, ignore_index=True).to_csv(args.output_dir / "direct_positive_gene_ranks.tsv", sep="\t", index=False)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
