#!/usr/bin/env python3
"""Harmonize GEO CLASH and CLEAR-CLIP interactions to miRNA--gene pairs."""

from __future__ import annotations

import argparse
import gzip
import re
from pathlib import Path

import pandas as pd


MIRNA_ALIASES = {
    "miR-16": "hsa-miR-16-5p",
    "miR-17": "hsa-miR-17-5p",
    "miR-25": "hsa-miR-25-3p",
    "miR-92a": "hsa-miR-92a-3p",
    "miR-124": "hsa-miR-124-3p",
    "miR-155": "hsa-miR-155-5p",
    "miR-186": "hsa-miR-186-5p",
}


def load_ensgene(path: Path):
    models = {}
    with gzip.open(path, "rt") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 13:
                continue
            transcript, strand = fields[1], fields[3]
            cds_start, cds_end = int(fields[6]), int(fields[7])
            starts = [int(x) for x in fields[9].rstrip(",").split(",")]
            ends = [int(x) for x in fields[10].rstrip(",").split(",")]
            if cds_start == cds_end:
                models[transcript] = ("noncoding_exonic", None, None)
                continue
            if strand == "+":
                five_end = sum(max(0, min(end, cds_start) - start) for start, end in zip(starts, ends) if start < cds_start)
                coding_end = sum(max(0, min(end, cds_end) - start) for start, end in zip(starts, ends) if start < cds_end)
            else:
                five_end = sum(max(0, end - max(start, cds_end)) for start, end in zip(starts, ends) if end > cds_end)
                coding_end = sum(max(0, end - max(start, cds_start)) for start, end in zip(starts, ends) if end > cds_start)
            models[transcript] = ("coding", five_end, coding_end)
    return models


def transcript_region(model, start: int, end: int) -> str:
    if model is None:
        return "mRNA_region_unmapped"
    kind, five_end, coding_end = model
    if kind == "noncoding_exonic":
        return kind
    midpoint = (start + end) / 2
    if midpoint <= five_end:
        return "5'UTR"
    if midpoint > coding_end:
        return "3'UTR"
    return "CDS"


def parse_clash(raw_dir: Path, ensgene: Path) -> pd.DataFrame:
    models = load_ensgene(ensgene)
    records = []
    for path in sorted(raw_dir.glob("*.hyb.txt.gz")):
        with gzip.open(path, "rt", errors="replace") as handle:
            for line in handle:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 15:
                    continue
                left, right = fields[3], fields[9]
                mir_field = left if "MirBase_" in left else right if "MirBase_" in right else None
                gene_field = right if "_mRNA" in right else left if "_mRNA" in left else None
                if not mir_field or not gene_field:
                    continue
                mir_match = re.search(r"MirBase_(miR-[^_]+)_microRNA", mir_field)
                gene_match = re.match(r"(ENSG\d+)_(ENST\d+)_([^_]+)_mRNA", gene_field)
                if not mir_match or not gene_match:
                    continue
                source_mirna = mir_match.group(1)
                current_mirna = MIRNA_ALIASES.get(source_mirna)
                if current_mirna is None:
                    continue
                count_match = re.search(r"_(\d+)$", fields[0])
                gene_on_right = gene_field == right
                transcript_start = int(fields[12] if gene_on_right else fields[6])
                transcript_end = int(fields[13] if gene_on_right else fields[7])
                transcript = gene_match.group(2)
                records.append({
                    "dataset": "GSE50452_CLASH",
                    "cell_line": "HEK293",
                    "mirna": current_mirna,
                    "source_mirna": source_mirna,
                    "gene": gene_match.group(3),
                    "ensembl_gene": gene_match.group(1),
                    "ensembl_transcript": transcript,
                    "region": transcript_region(models.get(transcript), transcript_start, transcript_end),
                    "transcript_start": transcript_start,
                    "transcript_end": transcript_end,
                    "support_count": int(count_match.group(1)) if count_match else 1,
                    "source_record": fields[0],
                    "source_file": path.name,
                })
    return pd.DataFrame.from_records(records)


def parse_clear_clip(workbook: Path) -> pd.DataFrame:
    frame = pd.read_excel(workbook, sheet_name="Huh-7.5 chimeras", header=23)
    frame = frame[frame["miRNA"].isin(set(MIRNA_ALIASES.values())) & frame["gene.symbol"].notna()].copy()
    return pd.DataFrame({
        "dataset": "GSE73057_CLEAR_CLIP",
        "cell_line": "Huh7.5",
        "mirna": frame["miRNA"].astype(str),
        "source_mirna": frame["miRNA"].astype(str),
        "gene": frame["gene.symbol"].astype(str),
        "ensembl_gene": "",
        "ensembl_transcript": "",
        "region": frame["region"].astype(str),
        "transcript_start": "",
        "transcript_end": "",
        "support_count": pd.to_numeric(frame["N"], errors="coerce").fillna(1).astype(int),
        "source_record": frame["cluster.ID"].astype(str),
        "source_file": workbook.name,
    })


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--clash-dir", required=True, type=Path)
    parser.add_argument("--clear-workbook", required=True, type=Path)
    parser.add_argument("--ensgene", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    combined = pd.concat([parse_clash(args.clash_dir, args.ensgene), parse_clear_clip(args.clear_workbook)], ignore_index=True)
    combined = combined.sort_values(["dataset", "mirna", "gene", "source_record"])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    combined.to_csv(args.output, sep="\t", index=False)
    print(combined.groupby(["dataset", "mirna"]).agg(records=("gene", "size"), genes=("gene", "nunique")).to_string())


if __name__ == "__main__":
    main()
