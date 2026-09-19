#!/usr/bin/env python3
"""Convert GSE52530 processed expression matrices to analysis-ready long TSV.

The GEO files have no header. Their fields are RefSeq accession, transcript
length, mock replicate values, miR-124 replicate values, and miR-155 replicate
values. Replicate values inside each condition are comma-separated. Keeping the
data long avoids silently choosing a pseudocount or collapsing replicates.
"""

from __future__ import annotations

import argparse
import csv
import gzip
from pathlib import Path
from typing import TextIO


CONDITION_LABELS = {
    "Puc19": "mock",
    "Puc19_124": "mock_miR-124",
    "Puc19_155": "mock_miR-155",
    "miR-124": "miR-124",
    "miR-155": "miR-155",
}


def _open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


def convert(input_path: Path, output_path: Path, cell_line: str) -> int:
    """Write long-format records and return the number of input transcripts."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    transcript_count = 0
    with _open_text(input_path) as source, output_path.open(
        "w", encoding="utf-8", newline=""
    ) as destination:
        reader = csv.reader(source, delimiter="\t")
        writer = csv.writer(destination, delimiter="\t", lineterminator="\n")
        writer.writerow([
            "cell_line",
            "refseq_id",
            "transcript_length",
            "condition",
            "replicate",
            "sample_id",
            "rpkm",
        ])

        try:
            header = next(reader)
        except StopIteration as exc:
            raise ValueError(f"{input_path}: empty input") from exc
        if len(header) < 3 or header[:2] != ["#geneID", "geneLen"]:
            raise ValueError(
                f"{input_path}: expected #geneID/geneLen header; use the quantile-normalized file"
            )
        groups = []
        for group_text in header[2:]:
            sample_ids = group_text.split(",")
            group_key = sample_ids[0].rsplit("_", 1)[0]
            if group_key not in CONDITION_LABELS:
                raise ValueError(f"{input_path}: unsupported sample group {group_key!r}")
            groups.append((CONDITION_LABELS[group_key], sample_ids))

        for line_number, row in enumerate(reader, start=2):
            if not row:
                continue
            if len(row) != len(header):
                raise ValueError(
                    f"{input_path}:{line_number}: expected {len(header)} tab-separated fields, "
                    f"found {len(row)}"
                )
            refseq_id, length_text, *condition_fields = row
            try:
                transcript_length = int(length_text)
            except ValueError as exc:
                raise ValueError(
                    f"{input_path}:{line_number}: invalid transcript length {length_text!r}"
                ) from exc

            for (condition, sample_ids), values_text in zip(
                groups, condition_fields, strict=True
            ):
                values = values_text.split(",")
                if len(values) != len(sample_ids) or any(value == "" for value in values):
                    raise ValueError(
                        f"{input_path}:{line_number}: sample/value count mismatch for {condition}"
                    )
                for replicate, (sample_id, value_text) in enumerate(
                    zip(sample_ids, values, strict=True), start=1
                ):
                    try:
                        rpkm = float(value_text)
                    except ValueError as exc:
                        raise ValueError(
                            f"{input_path}:{line_number}: invalid RPKM {value_text!r}"
                        ) from exc
                    writer.writerow(
                        [
                            cell_line,
                            refseq_id,
                            transcript_length,
                            condition,
                            replicate,
                            sample_id,
                            rpkm,
                        ]
                    )
            transcript_count += 1
    return transcript_count


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="GSE52530 *.expData[.qn].txt.gz file")
    parser.add_argument("output", type=Path, help="Output long-format TSV")
    parser.add_argument("--cell-line", required=True, help="Cell-line label stored in output")
    args = parser.parse_args()
    count = convert(args.input, args.output, args.cell_line)
    print(f"Wrote {count} transcripts to {args.output}")


if __name__ == "__main__":
    main()
