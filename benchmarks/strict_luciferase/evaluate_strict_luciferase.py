#!/usr/bin/env python3
"""Join manually audited WT/mutant 3'UTR reporters to frozen APAmiRank ranks.

This deliberately reports descriptive pilot summaries only. The audited pairs were
assembled after the existing reporter benchmark was inspected, so they cannot be
used as an unbiased estimate of generalization performance.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


STRICT_REQUIRED = {
    "target_region": "3UTR",
    "wt_repressed": "yes",
    "mutant_rescue": "yes",
    "mirna_specificity": "exact",
}


def validate_audit(audit: pd.DataFrame) -> None:
    required = {
        "miRNA", "gene", "target_region", "wt_repressed", "mutant_rescue",
        "mirna_specificity", "verification_level", "eligibility", "primary_url",
    }
    missing = sorted(required - set(audit.columns))
    if missing:
        raise ValueError(f"audit table is missing required columns: {missing}")
    if audit[["miRNA", "gene"]].isna().any().any():
        raise ValueError("miRNA and gene must be present for every audit row")
    strict = audit.eligibility.isin(["strict", "corroborated_strict"])
    for column, expected in STRICT_REQUIRED.items():
        bad = audit.loc[strict & audit[column].ne(expected), ["miRNA", "gene", column]]
        if not bad.empty:
            raise ValueError(f"eligible rows violate {column}={expected}: {bad.to_dict('records')}")
    if audit.loc[strict, ["miRNA", "gene"]].duplicated().any():
        raise ValueError("eligible miRNA-gene pairs must be unique")


def select_scope(audit: pd.DataFrame, scope: str) -> pd.DataFrame:
    if scope == "full_text_strict":
        keep = audit.eligibility.eq("strict") & audit.verification_level.eq("full_text_verified")
    elif scope == "plus_construct_corroborated":
        keep = audit.eligibility.isin(["strict", "corroborated_strict"])
    else:
        raise ValueError(f"unknown scope: {scope}")
    return audit.loc[keep].copy()


def summarize(joined: pd.DataFrame, source_n: int, scope: str) -> dict:
    observed = joined.consensus_apa_exposure_percentile.dropna()
    delta = (
        joined.consensus_apa_exposure_percentile
        - joined.consensus_unweighted_additive_percentile
    ).dropna()
    return {
        "scope": scope,
        "source_pairs": source_n,
        "observable_pairs": int(len(observed)),
        "candidate_coverage": len(observed) / source_n if source_n else np.nan,
        "mean_apa_percentile": observed.mean() if len(observed) else np.nan,
        "median_apa_percentile": observed.median() if len(observed) else np.nan,
        "top_decile_pairs": int((observed >= 0.90).sum()),
        "top_decile_fraction": (observed >= 0.90).mean() if len(observed) else np.nan,
        "top_quintile_pairs": int((observed >= 0.80).sum()),
        "top_quintile_fraction": (observed >= 0.80).mean() if len(observed) else np.nan,
        "median_apa_minus_unweighted_percentile": delta.median() if len(delta) else np.nan,
        "inferential_test": "not_run_post_selection_small_n",
    }


def evaluate(audit: pd.DataFrame, ranks: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    validate_audit(audit)
    scopes = ["full_text_strict", "plus_construct_corroborated"]
    details, summaries = [], []
    for scope in scopes:
        selected = select_scope(audit, scope)
        joined = selected.merge(ranks, on=["miRNA", "gene"], how="left", validate="one_to_one")
        joined.insert(0, "scope", scope)
        details.append(joined)
        summaries.append(summarize(joined, len(selected), scope))
    return pd.concat(details, ignore_index=True), pd.DataFrame(summaries)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--audit", required=True, type=Path)
    parser.add_argument("--ranks", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    audit = pd.read_csv(args.audit, sep="\t", dtype=str, keep_default_na=False)
    ranks = pd.read_csv(args.ranks, sep="\t")
    details, summary = evaluate(audit, ranks)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    details.to_csv(args.output_dir / "strict_luciferase_pair_ranks.tsv", sep="\t", index=False)
    summary.to_csv(args.output_dir / "strict_luciferase_summary.tsv", sep="\t", index=False)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
