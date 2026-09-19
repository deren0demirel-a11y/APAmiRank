#!/usr/bin/env python3
"""Initialize a rank-free paper-screening form from the locked manifest."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


SCREENING_FIELDS = [
    "screening_status", "article_title", "pmid", "pmcid", "doi", "primary_url",
    "article_status", "species", "mature_mirna_exact", "target_region",
    "reporter_vector", "cell_context", "wt_repression", "mutant_design",
    "mutant_rescue", "orthogonal_endogenous_support", "eligibility_decision",
    "exclusion_code", "evidence_note", "reviewer", "reviewed_on",
]


def initialize(manifest: pd.DataFrame) -> pd.DataFrame:
    form = manifest[["miRNA", "gene", "mti_ids", "selection_hash"]].copy()
    for field in SCREENING_FIELDS:
        form[field] = ""
    form["screening_status"] = "pending"
    form["eligibility_decision"] = "pending"
    return form


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    manifest = pd.read_csv(args.manifest, sep="\t", dtype=str)
    form = initialize(manifest)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    form.to_csv(args.output, sep="\t", index=False)
    print(f"initialized {len(form)} pending records")


if __name__ == "__main__":
    main()
