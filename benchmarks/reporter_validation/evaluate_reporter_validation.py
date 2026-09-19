#!/usr/bin/env python3
"""Evaluate consensus and cell-context APAmiRank ranks against reporter evidence."""

from __future__ import annotations

import argparse
from functools import reduce
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import hypergeom, mannwhitneyu


PANEL = [
    "hsa-miR-16-5p", "hsa-miR-17-5p", "hsa-miR-25-3p",
    "hsa-miR-92a-3p", "hsa-miR-124-3p", "hsa-miR-155-5p",
    "hsa-miR-186-5p",
]
CELLS = ["HeLa", "HEK293", "Huh7", "IMR90"]
EVIDENCE_SETS = {
    "primary_reporter_plus_wb_or_qpcr": "primary_reporter_plus_wb_or_qpcr",
    "reporter_any": None,
    "tier_a_or_b": "tier_a_or_b",
    "tier_a": "tier_a",
}


def percentile(rank: pd.Series, n: int) -> pd.Series:
    if n <= 1:
        return pd.Series(np.ones(len(rank)), index=rank.index)
    return 1.0 - (rank.astype(float) - 1.0) / (n - 1.0)


def load_cell_ranks(runs: Path, cell: str, mirna: str) -> pd.DataFrame:
    run = runs / f"{cell}__{mirna}"
    base = pd.read_csv(run / "03_gene_ranking.tsv", sep="\t")
    cond = pd.read_csv(run / "03c_condition_gene_ranking.tsv", sep="\t")
    sites = pd.read_csv(run / "02_site_ranking.tsv", sep="\t")
    n = len(base)
    result = base[["gene", "best_equal_weight_rank", "opportunity_adjusted_descriptive_rank"]].copy()
    result["best_site"] = percentile(result["best_equal_weight_rank"], n)
    result["opportunity_adjusted"] = percentile(result["opportunity_adjusted_descriptive_rank"], n)
    additive = sites.groupby("gene", as_index=False)["binding_evidence_score_equal_weight"].sum()
    additive["rank"] = additive["binding_evidence_score_equal_weight"].rank(method="min", ascending=False)
    additive["unweighted_additive"] = percentile(additive["rank"], n)
    cond = cond[["gene", "condition_gene_rank"]].copy()
    cond["apa_exposure"] = percentile(cond["condition_gene_rank"], n)
    return result[["gene", "best_site", "opportunity_adjusted"]].merge(
        additive[["gene", "unweighted_additive"]], on="gene", how="inner"
    ).merge(cond[["gene", "apa_exposure"]], on="gene", how="inner")


def bh_adjust(frame: pd.DataFrame, group_cols: list[str]) -> pd.DataFrame:
    frame = frame.copy()
    frame["mann_whitney_bh_q"] = np.nan
    for _, idx in frame.groupby(group_cols).groups.items():
        valid = frame.loc[idx, "mann_whitney_one_sided_p"].dropna().sort_values()
        if valid.empty:
            continue
        raw = valid.to_numpy() * len(valid) / np.arange(1, len(valid) + 1)
        adjusted = np.minimum.accumulate(raw[::-1])[::-1].clip(max=1.0)
        frame.loc[valid.index, "mann_whitney_bh_q"] = adjusted
    return frame


def metrics(ranks: pd.DataFrame, positive_genes: set[str], score_col: str) -> dict:
    n = len(ranks)
    positive = ranks.loc[ranks.gene.isin(positive_genes), score_col].dropna()
    unlabeled = ranks.loc[~ranks.gene.isin(positive_genes), score_col].dropna()
    if len(positive) and len(unlabeled):
        test = mannwhitneyu(positive, unlabeled, alternative="greater")
        auc = test.statistic / (len(positive) * len(unlabeled))
        p_value = test.pvalue
    else:
        auc = p_value = np.nan
    row = {
        "candidate_genes": n,
        "source_positive_genes": len(positive_genes),
        "overlapping_positive_genes": len(positive),
        "candidate_coverage_of_source_positives": len(positive) / len(positive_genes) if positive_genes else np.nan,
        "median_positive_percentile": positive.median() if len(positive) else np.nan,
        "positive_unlabeled_auc": auc,
        "mann_whitney_one_sided_p": p_value,
    }
    for fraction in (0.10, 0.20):
        top_n = max(1, int(np.ceil(fraction * n)))
        threshold = 1 - (top_n - 1) / (n - 1) if n > 1 else 1
        hits = int((positive >= threshold).sum())
        label = f"top_{int(fraction * 100)}"
        row[f"{label}_hits"] = hits
        row[f"{label}_source_recall"] = hits / len(positive_genes) if positive_genes else np.nan
        row[f"{label}_observable_recall"] = hits / len(positive) if len(positive) else np.nan
        row[f"{label}_enrichment"] = (hits / len(positive)) / (top_n / n) if len(positive) else np.nan
        row[f"{label}_hypergeom_p"] = hypergeom.sf(hits - 1, n, len(positive), top_n) if len(positive) else np.nan
    return row


def evidence_genes(evidence: pd.DataFrame, mirna: str, evidence_set: str) -> set[str]:
    subset = evidence[evidence.miRNA.eq(mirna)]
    flag = EVIDENCE_SETS[evidence_set]
    if flag is not None:
        subset = subset[subset[flag].astype(bool)]
    return set(subset.gene)


def auc_from_arrays(positive: np.ndarray, unlabeled: np.ndarray) -> float:
    ordered = np.sort(unlabeled)
    lower = np.searchsorted(ordered, positive, side="left")
    upper = np.searchsorted(ordered, positive, side="right")
    return float(np.mean((lower + 0.5 * (upper - lower)) / len(ordered)))


def bootstrap_auc_increment(ranks: pd.DataFrame, positive_genes: set[str], iterations: int, seed: int) -> dict:
    labels = ranks.gene.isin(positive_genes).to_numpy()
    pos = ranks.loc[labels, ["consensus_apa_exposure_percentile", "consensus_unweighted_additive_percentile"]].to_numpy()
    neg = ranks.loc[~labels, ["consensus_apa_exposure_percentile", "consensus_unweighted_additive_percentile"]].to_numpy()
    if len(pos) < 5 or len(neg) < 5:
        return {"apa_minus_unweighted_auc": np.nan, "bootstrap_ci_low": np.nan, "bootstrap_ci_high": np.nan, "bootstrap_one_sided_p": np.nan}
    observed = auc_from_arrays(pos[:, 0], neg[:, 0]) - auc_from_arrays(pos[:, 1], neg[:, 1])
    rng = np.random.default_rng(seed)
    draws = np.empty(iterations)
    for index in range(iterations):
        pos_sample = pos[rng.integers(0, len(pos), len(pos))]
        neg_sample = neg[rng.integers(0, len(neg), len(neg))]
        draws[index] = auc_from_arrays(pos_sample[:, 0], neg_sample[:, 0]) - auc_from_arrays(pos_sample[:, 1], neg_sample[:, 1])
    return {
        "apa_minus_unweighted_auc": observed,
        "bootstrap_ci_low": np.quantile(draws, 0.025),
        "bootstrap_ci_high": np.quantile(draws, 0.975),
        "bootstrap_one_sided_p": (1 + np.sum(draws <= 0)) / (iterations + 1),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evidence", required=True, type=Path)
    parser.add_argument("--runs", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    evidence = pd.read_csv(args.evidence, sep="\t")
    consensus_rows, cell_rows, details, increments = [], [], [], []
    for mirna in PANEL:
        by_cell = {cell: load_cell_ranks(args.runs, cell, mirna) for cell in CELLS}
        renamed = []
        for cell, table in by_cell.items():
            renamed.append(table.rename(columns={c: f"{c}__{cell}" for c in table.columns if c != "gene"}))
        outer = reduce(lambda left, right: left.merge(right, on="gene", how="outer"), renamed)
        outer["observed_contexts"] = outer[[f"best_site__{cell}" for cell in CELLS]].notna().sum(axis=1)
        universes = {
            "at_least_2_of_4_contexts": outer[outer.observed_contexts >= 2].copy(),
            "all_4_contexts": outer[outer.observed_contexts == 4].copy(),
        }
        for universe_scope, consensus in universes.items():
            for metric_name in ["best_site", "opportunity_adjusted", "unweighted_additive", "apa_exposure"]:
                cols = [f"{metric_name}__{cell}" for cell in CELLS]
                raw_col = f"consensus_{metric_name}_mean"
                rank_col = f"consensus_{metric_name}_percentile"
                consensus[raw_col] = consensus[cols].mean(axis=1, skipna=True)
                consensus[rank_col] = percentile(consensus[raw_col].rank(method="min", ascending=False), len(consensus))
            for evidence_set in EVIDENCE_SETS:
                positive_genes = evidence_genes(evidence, mirna, evidence_set)
                for metric_name in ["best_site", "opportunity_adjusted", "unweighted_additive", "apa_exposure"]:
                    row = metrics(consensus, positive_genes, f"consensus_{metric_name}_percentile")
                    row.update({"miRNA": mirna, "universe_scope": universe_scope, "evidence_set": evidence_set, "rank_metric": metric_name})
                    consensus_rows.append(row)
                increment = bootstrap_auc_increment(consensus, positive_genes, iterations=5000, seed=20260908)
                increment.update({
                    "miRNA": mirna,
                    "universe_scope": universe_scope,
                    "evidence_set": evidence_set,
                    "overlapping_positive_genes": int(consensus.gene.isin(positive_genes).sum()),
                    "bootstrap_iterations": 5000,
                })
                increments.append(increment)
                if evidence_set == "primary_reporter_plus_wb_or_qpcr" and universe_scope == "at_least_2_of_4_contexts":
                    positive_detail = consensus[consensus.gene.isin(positive_genes)].copy()
                    positive_detail.insert(0, "miRNA", mirna)
                    positive_detail.insert(1, "universe_scope", universe_scope)
                    details.append(positive_detail)
        positive_genes = evidence_genes(evidence, mirna, "primary_reporter_plus_wb_or_qpcr")
        for cell, table in by_cell.items():
            row = metrics(table, positive_genes, "apa_exposure")
            row.update({"miRNA": mirna, "cell_line": cell, "evidence_set": "primary_reporter_plus_wb_or_qpcr", "rank_metric": "apa_exposure"})
            cell_rows.append(row)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    consensus_summary = bh_adjust(pd.DataFrame(consensus_rows), ["universe_scope", "evidence_set", "rank_metric"])
    cell_summary = bh_adjust(pd.DataFrame(cell_rows), ["evidence_set", "rank_metric"])
    increment_summary = pd.DataFrame(increments)
    increment_summary["bootstrap_bh_q"] = np.nan
    for _, idx in increment_summary.groupby(["universe_scope", "evidence_set"]).groups.items():
        valid = increment_summary.loc[idx, "bootstrap_one_sided_p"].dropna().sort_values()
        if valid.empty:
            continue
        raw = valid.to_numpy() * len(valid) / np.arange(1, len(valid) + 1)
        adjusted = np.minimum.accumulate(raw[::-1])[::-1].clip(max=1.0)
        increment_summary.loc[valid.index, "bootstrap_bh_q"] = adjusted
    consensus_summary.to_csv(args.output_dir / "reporter_consensus_summary.tsv", sep="\t", index=False)
    cell_summary.to_csv(args.output_dir / "reporter_cell_context_summary.tsv", sep="\t", index=False)
    increment_summary.to_csv(args.output_dir / "apa_increment_bootstrap.tsv", sep="\t", index=False)
    pd.concat(details, ignore_index=True).to_csv(args.output_dir / "primary_reporter_positive_consensus_ranks.tsv", sep="\t", index=False)
    primary = consensus_summary[
        (consensus_summary.evidence_set == "primary_reporter_plus_wb_or_qpcr")
        & (consensus_summary.universe_scope == "at_least_2_of_4_contexts")
        & (consensus_summary.rank_metric == "apa_exposure")
    ]
    print(primary.to_string(index=False))


if __name__ == "__main__":
    main()
