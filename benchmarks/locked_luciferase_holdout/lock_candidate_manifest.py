#!/usr/bin/env python3
"""Create a deterministic, rank-independent luciferase screening manifest."""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

import pandas as pd


SEED = "APAmiRank-locked-luciferase-holdout-v1-2026-09-08"
PANEL = [
    "hsa-miR-16-5p", "hsa-miR-17-5p", "hsa-miR-25-3p",
    "hsa-miR-92a-3p", "hsa-miR-124-3p", "hsa-miR-155-5p",
    "hsa-miR-186-5p",
]


def stable_hash(mirna: str, gene: str) -> str:
    payload = f"{SEED}\t{mirna}\t{gene}".encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def select_candidates(
    evidence: pd.DataFrame,
    prior_audit: pd.DataFrame,
    per_mirna: int = 5,
) -> pd.DataFrame:
    required = {
        "miRNA", "gene", "mti_ids", "reporter_assay", "western_blot", "qpcr",
        "tier_a_or_b", "max_papers", "primary_reporter_plus_wb_or_qpcr",
    }
    missing = sorted(required - set(evidence.columns))
    if missing:
        raise ValueError(f"evidence table is missing required columns: {missing}")

    excluded_pairs = set(map(tuple, prior_audit[["miRNA", "gene"]].to_numpy()))
    eligible = evidence[
        evidence.miRNA.isin(PANEL)
        & evidence.reporter_assay.astype(bool)
        & evidence.primary_reporter_plus_wb_or_qpcr.astype(bool)
        & evidence.tier_a_or_b.astype(bool)
        & (pd.to_numeric(evidence.max_papers, errors="coerce") >= 1)
    ].copy()
    eligible = eligible[
        ~eligible[["miRNA", "gene"]].apply(tuple, axis=1).isin(excluded_pairs)
    ].copy()
    eligible["selection_hash"] = [
        stable_hash(mirna, gene) for mirna, gene in eligible[["miRNA", "gene"]].to_numpy()
    ]
    eligible = eligible.sort_values(["miRNA", "selection_hash", "gene"])
    selected = eligible.groupby("miRNA", sort=False, group_keys=False).head(per_mirna).copy()
    selected["eligible_pool_size"] = selected.miRNA.map(eligible.groupby("miRNA").size())
    selected["requested_per_mirna"] = per_mirna
    selected["selection_seed"] = SEED
    selected["selection_rule"] = (
        "human panel; reporter; WB_or_qPCR; Tier_A_or_B; max_papers>=1; "
        "exclude_prior_audit; lowest_SHA256"
    )
    columns = [
        "miRNA", "gene", "mti_ids", "reporter_assay", "western_blot", "qpcr",
        "tier_a_or_b", "max_papers", "eligible_pool_size", "requested_per_mirna",
        "selection_hash", "selection_seed", "selection_rule",
    ]
    return selected[columns].reset_index(drop=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evidence", required=True, type=Path)
    parser.add_argument("--prior-audit", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--per-mirna", type=int, default=5)
    args = parser.parse_args()
    evidence = pd.read_csv(args.evidence, sep="\t")
    prior = pd.read_csv(args.prior_audit, sep="\t", dtype=str, keep_default_na=False)
    selected = select_candidates(evidence, prior, args.per_mirna)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    selected.to_csv(args.output, sep="\t", index=False)
    print(selected.groupby("miRNA").size().to_string())


if __name__ == "__main__":
    main()
