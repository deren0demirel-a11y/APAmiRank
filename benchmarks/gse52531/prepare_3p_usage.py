#!/usr/bin/env python3
"""Map GSE52527 hg19 3P-seq clusters to conservative RefSeq tandem 3'UTRs.

The GEO BED records span mapped 3P-tag clusters. The transcriptional 3' boundary
is used as the PAS coordinate (BED end for '+' and BED start for '-'). Only PASs
inside an annotated, contiguous terminal-exon 3'UTR are retained. Ambiguous
gene/locus assignments are excluded rather than resolved heuristically.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import TextIO


@dataclass(frozen=True)
class Locus:
    gene: str
    chrom: str
    strand: str
    anchor: int
    distal: int
    accessions: tuple[str, ...]

    @property
    def start0(self) -> int:
        return self.anchor if self.strand == "+" else self.distal

    @property
    def end(self) -> int:
        return self.distal if self.strand == "+" else self.anchor

    @property
    def key(self) -> tuple[str, str, str, int]:
        return (self.gene, self.chrom, self.strand, self.anchor)


def _open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


def load_refgene(path: Path) -> list[Locus]:
    grouped: dict[tuple[str, str, str, int], dict[str, object]] = {}
    with _open_text(path) as handle:
        reader = csv.reader(handle, delimiter="\t")
        for line_number, row in enumerate(reader, start=1):
            if len(row) < 13:
                raise ValueError(f"{path}:{line_number}: expected at least 13 fields")
            accession, chrom, strand = row[1], row[2], row[3]
            tx_start, tx_end = int(row[4]), int(row[5])
            cds_start, cds_end = int(row[6]), int(row[7])
            exon_starts = [int(value) for value in row[9].rstrip(",").split(",")]
            exon_ends = [int(value) for value in row[10].rstrip(",").split(",")]
            gene = row[12]
            if strand not in {"+", "-"} or cds_start >= cds_end:
                continue
            if not chrom.startswith("chr") or "_" in chrom:
                continue

            if strand == "+":
                terminal_start, terminal_end = exon_starts[-1], exon_ends[-1]
                if not (terminal_start <= cds_end < terminal_end):
                    continue
                anchor, distal = cds_end, tx_end
            else:
                terminal_start, terminal_end = exon_starts[0], exon_ends[0]
                if not (terminal_start < cds_start <= terminal_end):
                    continue
                anchor, distal = cds_start, tx_start
            if anchor == distal:
                continue

            key = (gene, chrom, strand, anchor)
            if key not in grouped:
                grouped[key] = {"distal": distal, "accessions": set()}
            else:
                current = int(grouped[key]["distal"])
                is_more_distal = distal > current if strand == "+" else distal < current
                if is_more_distal:
                    grouped[key]["distal"] = distal
                    grouped[key]["accessions"] = set()
            if distal == int(grouped[key]["distal"]):
                grouped[key]["accessions"].add(accession)

    return [
        Locus(
            gene=key[0],
            chrom=key[1],
            strand=key[2],
            anchor=key[3],
            distal=int(values["distal"]),
            accessions=tuple(sorted(values["accessions"])),
        )
        for key, values in sorted(grouped.items())
    ]


def build_bin_index(loci: list[Locus], bin_size: int = 100_000):
    index: dict[tuple[str, str, int], list[Locus]] = defaultdict(list)
    for locus in loci:
        for bin_id in range(locus.start0 // bin_size, locus.end // bin_size + 1):
            index[(locus.chrom, locus.strand, bin_id)].append(locus)
    return index


def map_clusters(
    bed_path: Path,
    loci: list[Locus],
    min_count: float = 2.0,
    bin_size: int = 100_000,
):
    index = build_bin_index(loci, bin_size)
    counts: dict[tuple[tuple[str, str, str, int], int], float] = defaultdict(float)
    qc = defaultdict(int)
    locus_lookup = {locus.key: locus for locus in loci}

    with _open_text(bed_path) as handle:
        reader = csv.reader(handle, delimiter="\t")
        for line_number, row in enumerate(reader, start=1):
            qc["input_clusters"] += 1
            if len(row) < 6:
                raise ValueError(f"{bed_path}:{line_number}: expected at least 6 BED fields")
            chrom, strand = row[0], row[5]
            start0, end = int(row[1]), int(row[2])
            try:
                count = float(row[3])
            except ValueError as exc:
                raise ValueError(f"{bed_path}:{line_number}: invalid cluster count") from exc
            if not math.isfinite(count) or count < min_count:
                qc["below_min_count"] += 1
                continue
            if strand not in {"+", "-"} or start0 >= end:
                qc["invalid_bed_record"] += 1
                continue
            pas = end if strand == "+" else start0
            candidates = []
            for locus in index.get((chrom, strand, pas // bin_size), []):
                if locus.strand == "+" and locus.anchor < pas <= locus.distal:
                    candidates.append(locus)
                elif locus.strand == "-" and locus.distal <= pas < locus.anchor:
                    candidates.append(locus)
            unique = {candidate.key: candidate for candidate in candidates}
            if not unique:
                qc["outside_conservative_terminal_utr"] += 1
                continue
            if len(unique) != 1:
                qc["ambiguous_locus"] += 1
                continue
            locus = next(iter(unique.values()))
            counts[(locus.key, pas)] += count
            qc["mapped_cluster_records"] += 1

    return counts, locus_lookup, qc


def write_outputs(
    counts,
    locus_lookup,
    condition: str,
    bed_out: Path,
    usage_out: Path,
    map_out: Path,
    qc_out: Path,
    qc,
    min_pas_per_locus: int = 2,
    min_locus_count: float = 10.0,
    min_isoform_fraction: float = 0.01,
):
    by_locus = defaultdict(list)
    for (locus_key, pas), count in counts.items():
        by_locus[locus_key].append((pas, count))

    retained = {}
    low_depth_loci = 0
    non_apa_loci = 0
    for key, values in by_locus.items():
        total = sum(count for _, count in values)
        if total < min_locus_count:
            low_depth_loci += 1
            continue
        supported_pas = sum(count / total >= min_isoform_fraction for _, count in values)
        if supported_pas < min_pas_per_locus:
            non_apa_loci += 1
            continue
        # Retain minor PASs in qualifying loci so their tag mass remains in the
        # denominator and usage fractions sum to one without renormalization.
        retained[key] = values
    bed_out.parent.mkdir(parents=True, exist_ok=True)
    usage_out.parent.mkdir(parents=True, exist_ok=True)
    map_out.parent.mkdir(parents=True, exist_ok=True)

    with bed_out.open("w", encoding="utf-8", newline="") as bed_handle, usage_out.open(
        "w", encoding="utf-8", newline=""
    ) as usage_handle, map_out.open("w", encoding="utf-8", newline="") as map_handle:
        bed_writer = csv.writer(bed_handle, delimiter="\t", lineterminator="\n")
        usage_writer = csv.writer(usage_handle, delimiter="\t", lineterminator="\n")
        map_writer = csv.writer(map_handle, delimiter="\t", lineterminator="\n")
        usage_writer.writerow(["isoform_id", "condition", "usage_fraction"])
        map_writer.writerow(
            ["isoform_id", "gene", "chrom", "strand", "anchor", "pas", "cluster_count", "refseq_accessions"]
        )

        for locus_key in sorted(retained):
            locus = locus_lookup[locus_key]
            sites = sorted(
                retained[locus_key], key=lambda item: item[0], reverse=locus.strand == "-"
            )
            total = sum(count for _, count in sites)
            for pas, count in sites:
                strand_label = "plus" if locus.strand == "+" else "minus"
                isoform_id = (
                    f"{locus.gene}__{locus.chrom}_{strand_label}_{locus.anchor}__PAS{pas}"
                )
                start0, end = (
                    (locus.anchor, pas) if locus.strand == "+" else (pas, locus.anchor)
                )
                bed_writer.writerow(
                    [locus.chrom, start0, end, isoform_id, f"{count:g}", locus.strand, locus.gene]
                )
                usage_writer.writerow([isoform_id, condition, f"{count / total:.12g}"])
                map_writer.writerow(
                    [
                        isoform_id,
                        locus.gene,
                        locus.chrom,
                        locus.strand,
                        locus.anchor,
                        pas,
                        f"{count:g}",
                        ",".join(locus.accessions),
                    ]
                )

    qc = dict(qc)
    qc.update(
        {
            "mapped_unique_pas": len(counts),
            "mapped_loci_before_apa_filter": len(by_locus),
            "retained_apa_loci": len(retained),
            "retained_isoforms": sum(len(values) for values in retained.values()),
            "min_pas_per_locus": min_pas_per_locus,
            "min_locus_count": min_locus_count,
            "min_isoform_fraction_for_apa": min_isoform_fraction,
            "low_depth_loci": low_depth_loci,
            "loci_without_two_supported_pas": non_apa_loci,
            "condition": condition,
            "pas_definition": "BED end for plus; BED start for minus",
            "assignment_policy": "single conservative RefSeq terminal-exon 3UTR locus only",
        }
    )
    qc_out.parent.mkdir(parents=True, exist_ok=True)
    qc_out.write_text(json.dumps(qc, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return qc


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bed", required=True, type=Path)
    parser.add_argument("--refgene", required=True, type=Path)
    parser.add_argument("--condition", required=True)
    parser.add_argument("--bed-out", required=True, type=Path)
    parser.add_argument("--usage-out", required=True, type=Path)
    parser.add_argument("--map-out", required=True, type=Path)
    parser.add_argument("--qc-out", required=True, type=Path)
    parser.add_argument("--min-count", type=float, default=2.0)
    parser.add_argument("--min-pas-per-locus", type=int, default=2)
    parser.add_argument("--min-locus-count", type=float, default=10.0)
    parser.add_argument("--min-isoform-fraction", type=float, default=0.01)
    args = parser.parse_args()

    loci = load_refgene(args.refgene)
    counts, locus_lookup, qc = map_clusters(args.bed, loci, args.min_count)
    report = write_outputs(
        counts,
        locus_lookup,
        args.condition,
        args.bed_out,
        args.usage_out,
        args.map_out,
        args.qc_out,
        qc,
        args.min_pas_per_locus,
        args.min_locus_count,
        args.min_isoform_fraction,
    )
    print(json.dumps(report, sort_keys=True))


if __name__ == "__main__":
    main()
