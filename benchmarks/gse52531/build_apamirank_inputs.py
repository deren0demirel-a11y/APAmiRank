#!/usr/bin/env python3
"""Build APAmiRank metadata and strand-correct FASTA from benchmark BED7."""

from __future__ import annotations

import argparse
import csv
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path


def parse_bed7(path: Path):
    rows = []
    seen = set()
    with path.open(encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for line_number, row in enumerate(reader, start=1):
            if len(row) < 7:
                raise ValueError(f"{path}:{line_number}: expected benchmark BED7")
            chrom, start_text, end_text, isoform_id, score, strand, gene = row[:7]
            start0, end = int(start_text), int(end_text)
            if isoform_id in seen:
                raise ValueError(f"{path}:{line_number}: duplicate isoform ID {isoform_id}")
            if strand not in {"+", "-"} or start0 >= end:
                raise ValueError(f"{path}:{line_number}: invalid interval or strand")
            seen.add(isoform_id)
            anchor = start0 if strand == "+" else end
            rows.append(
                {
                    "chrom": chrom,
                    "start0": start0,
                    "end": end,
                    "isoform_id": isoform_id,
                    "score": score,
                    "strand": strand,
                    "gene": gene,
                    "anchor": anchor,
                }
            )
    return rows


def write_metadata(rows, output: Path):
    groups = defaultdict(list)
    for row in rows:
        groups[(row["gene"], row["chrom"], row["strand"], row["anchor"])].append(row)
    metadata = []
    for locus_number, (key, values) in enumerate(sorted(groups.items()), start=1):
        gene, _chrom, _strand, _anchor = key
        locus_id = f"{gene}_L{locus_number}"
        values.sort(key=lambda row: (row["end"] - row["start0"], row["isoform_id"]))
        for rank, row in enumerate(values, start=1):
            metadata.append(
                {
                    "isoform_id": row["isoform_id"],
                    "gene": row["gene"],
                    "locus_id": locus_id,
                    "pas_rank": rank,
                    "chrom": row["chrom"],
                    "start0": row["start0"],
                    "end": row["end"],
                    "strand": row["strand"],
                }
            )
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(metadata[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(metadata)
    return metadata


def read_fasta_lengths(path: Path):
    lengths = {}
    current = None
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                current = line[1:].split()[0]
                if current in lengths:
                    raise ValueError(f"duplicate FASTA ID {current}")
                lengths[current] = 0
            elif line:
                if current is None:
                    raise ValueError("FASTA sequence before first header")
                lengths[current] += len(line)
    return lengths


def build(bed7: Path, twobit: Path, twobittofa: Path, metadata_out: Path, fasta_out: Path):
    rows = parse_bed7(bed7)
    if not rows:
        raise ValueError(f"{bed7}: no records")
    metadata = write_metadata(rows, metadata_out)
    fasta_out.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=fasta_out.parent) as temp_dir:
        bed6 = Path(temp_dir) / "intervals.bed"
        temporary_fasta = Path(temp_dir) / "isoforms.fa"
        with bed6.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            for row in rows:
                writer.writerow(
                    [
                        row["chrom"],
                        row["start0"],
                        row["end"],
                        row["isoform_id"],
                        row["score"],
                        row["strand"],
                    ]
                )
        subprocess.run(
            [str(twobittofa), str(twobit), f"-bed={bed6}", str(temporary_fasta)],
            check=True,
        )
        lengths = read_fasta_lengths(temporary_fasta)
        expected = {row["isoform_id"]: row["end"] - row["start0"] for row in rows}
        if lengths != expected:
            missing = sorted(set(expected) - set(lengths))[:5]
            mismatched = sorted(
                key for key in expected.keys() & lengths.keys() if expected[key] != lengths[key]
            )[:5]
            raise ValueError(
                f"sequence extraction mismatch; missing={missing}, length_mismatch={mismatched}"
            )
        temporary_fasta.replace(fasta_out)
    return metadata


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bed7", required=True, type=Path)
    parser.add_argument("--twobit", required=True, type=Path)
    parser.add_argument("--twobittofa", required=True, type=Path)
    parser.add_argument("--metadata-out", required=True, type=Path)
    parser.add_argument("--fasta-out", required=True, type=Path)
    args = parser.parse_args()
    metadata = build(
        args.bed7, args.twobit, args.twobittofa, args.metadata_out, args.fasta_out
    )
    print(f"Wrote {len(metadata)} isoforms to {args.metadata_out} and {args.fasta_out}")


if __name__ == "__main__":
    main()
