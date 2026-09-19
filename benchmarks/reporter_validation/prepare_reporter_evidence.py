#!/usr/bin/env python3
"""Harmonize a dated human miRTarBase Reporter Assay export for APAmiRank."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


PANEL = {
    "hsa-miR-16-5p", "hsa-miR-17-5p", "hsa-miR-25-3p",
    "hsa-miR-92a-3p", "hsa-miR-124-3p", "hsa-miR-155-5p",
    "hsa-miR-186-5p",
}
TIER_ORDER = {"Tier A: Strong": 3, "Tier B: High": 2, "Tier C: Middle": 1, "Tier D: Low": 0}


def join_unique(values):
    return ";".join(sorted({str(x) for x in values if pd.notna(x)}))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    data = pd.read_csv(args.input, encoding="utf-8-sig")
    data = data[
        data["Species (miRNA)"].eq("Homo sapiens")
        & data["Species (Target)"].eq("Homo sapiens")
        & data["miRNA"].isin(PANEL)
        & data["Reporter Assay"].eq(1)
    ].copy()
    data["gene"] = data["Target Gene"].astype(str).str.strip().str.upper()
    data["tier_value"] = data["Score"].map(TIER_ORDER).fillna(-1).astype(int)
    data["orthogonal_low_throughput"] = ((data["Western Blot"] == 1) | (data["qPCR"] == 1)).astype(int)
    grouped = data.groupby(["miRNA", "gene"], as_index=False).agg(
        mti_ids=("MTI ID", join_unique),
        reporter_assay=("Reporter Assay", "max"),
        western_blot=("Western Blot", "max"),
        qpcr=("qPCR", "max"),
        clip_seq=("CLIP-Seq", "max"),
        clash=("CLASH", "max"),
        orthogonal_low_throughput=("orthogonal_low_throughput", "max"),
        best_tier_value=("tier_value", "max"),
        max_papers=("Papers", "max"),
    )
    grouped["evidence_set"] = "reporter_any"
    grouped["primary_reporter_plus_wb_or_qpcr"] = grouped["orthogonal_low_throughput"].astype(bool)
    grouped["tier_a_or_b"] = grouped["best_tier_value"] >= 2
    grouped["tier_a"] = grouped["best_tier_value"] >= 3
    args.output.parent.mkdir(parents=True, exist_ok=True)
    grouped.sort_values(["miRNA", "gene"]).to_csv(args.output, sep="\t", index=False)
    print(grouped.groupby("miRNA").agg(
        reporter_genes=("gene", "size"),
        primary_genes=("primary_reporter_plus_wb_or_qpcr", "sum"),
        tier_a_genes=("tier_a", "sum"),
    ).to_string())


if __name__ == "__main__":
    main()
