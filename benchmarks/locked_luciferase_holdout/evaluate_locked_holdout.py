#!/usr/bin/env python3
"""Fail-closed evaluation of a completed, checksummed luciferase holdout."""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

import numpy as np
import pandas as pd


INCLUDED_REQUIRED = [
    "article_title", "primary_url", "species", "mature_mirna_exact",
    "target_region", "reporter_vector", "wt_repression", "mutant_design",
    "mutant_rescue", "evidence_note", "reviewer", "reviewed_on",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def validate_lock(manifest: pd.DataFrame, screening: pd.DataFrame, screening_path: Path, expected_sha256: str) -> None:
    if sha256(screening_path) != expected_sha256.lower():
        raise ValueError("screening form SHA-256 does not match the supplied lock")
    keys = ["miRNA", "gene", "selection_hash"]
    if screening[keys].duplicated().any():
        raise ValueError("screening form contains duplicate selected pairs")
    expected = set(map(tuple, manifest[keys].to_numpy()))
    observed = set(map(tuple, screening[keys].to_numpy()))
    if expected != observed:
        raise ValueError("screening form does not exactly match the locked manifest")
    if screening.screening_status.ne("complete").any() or screening.eligibility_decision.eq("pending").any():
        raise ValueError("screening is incomplete; rank join is blocked")
    valid_decisions = {"include", "exclude"}
    if not set(screening.eligibility_decision).issubset(valid_decisions):
        raise ValueError("eligibility_decision must be include or exclude")
    included = screening[screening.eligibility_decision.eq("include")]
    for column in INCLUDED_REQUIRED:
        if included[column].fillna("").str.strip().eq("").any():
            raise ValueError(f"included records require non-empty {column}")
    biological = {
        "species": "Homo sapiens", "mature_mirna_exact": "yes",
        "target_region": "3UTR", "wt_repression": "yes", "mutant_rescue": "yes",
    }
    for column, value in biological.items():
        if included[column].ne(value).any():
            raise ValueError(f"included records require {column}={value}")
    excluded = screening[screening.eligibility_decision.eq("exclude")]
    if excluded.exclusion_code.fillna("").str.strip().eq("").any():
        raise ValueError("excluded records require an exclusion_code")


def summarize(joined: pd.DataFrame, selected_n: int) -> pd.DataFrame:
    observed = joined.consensus_apa_exposure_percentile.dropna()
    delta = (
        joined.consensus_apa_exposure_percentile
        - joined.consensus_unweighted_additive_percentile
    ).dropna()
    row = {
        "selected_pairs": selected_n,
        "strict_included_pairs": len(joined),
        "observable_pairs": len(observed),
        "candidate_coverage": len(observed) / len(joined) if len(joined) else np.nan,
        "mean_apa_percentile": observed.mean() if len(observed) else np.nan,
        "median_apa_percentile": observed.median() if len(observed) else np.nan,
        "top_decile_fraction": (observed >= 0.90).mean() if len(observed) else np.nan,
        "top_quintile_fraction": (observed >= 0.80).mean() if len(observed) else np.nan,
        "median_apa_minus_unweighted_percentile": delta.median() if len(delta) else np.nan,
    }
    return pd.DataFrame([row])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--screening", required=True, type=Path)
    parser.add_argument("--screening-sha256", required=True)
    parser.add_argument("--ranks", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    manifest = pd.read_csv(args.manifest, sep="\t", dtype=str, keep_default_na=False)
    screening = pd.read_csv(args.screening, sep="\t", dtype=str, keep_default_na=False)
    validate_lock(manifest, screening, args.screening, args.screening_sha256)
    included = screening[screening.eligibility_decision.eq("include")].copy()
    ranks = pd.read_csv(args.ranks, sep="\t")
    joined = included.merge(ranks, on=["miRNA", "gene"], how="left", validate="one_to_one")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    joined.to_csv(args.output_dir / "locked_holdout_pair_ranks.tsv", sep="\t", index=False)
    summary = summarize(joined, len(manifest))
    summary.to_csv(args.output_dir / "locked_holdout_summary.tsv", sep="\t", index=False)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
