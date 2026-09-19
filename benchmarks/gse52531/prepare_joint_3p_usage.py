#!/usr/bin/env python3
"""Build a common hg19 PAS atlas and cell-line-specific usage matrix."""

from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path

from prepare_3p_usage import load_refgene, map_clusters


def cluster_positions(position_counts, max_span=30):
    """Complete-linkage-style genomic clustering with a bounded cluster span."""
    clusters = []
    for position in sorted(position_counts):
        if not clusters or position - clusters[-1][0] > max_span:
            clusters.append([position])
        else:
            clusters[-1].append(position)
    return clusters


def weighted_median(positions, weights):
    pairs = sorted(zip(positions, weights, strict=True))
    halfway = sum(weights) / 2
    cumulative = 0.0
    for position, weight in pairs:
        cumulative += weight
        if cumulative >= halfway:
            return position
    return pairs[-1][0]


def build_joint(
    samples,
    refgene,
    bed_out,
    usage_out,
    map_out,
    qc_out,
    min_cluster_count=2.0,
    min_condition_locus_count=10.0,
    min_isoform_fraction=0.01,
    max_pas_cluster_span=30,
):
    loci = load_refgene(refgene)
    lookup = {locus.key: locus for locus in loci}
    sample_counts = {}
    sample_qc = {}
    for condition, bed_path in samples:
        counts, _lookup, qc = map_clusters(bed_path, loci, min_cluster_count)
        sample_counts[condition] = counts
        sample_qc[condition] = dict(qc)

    positions_by_locus = defaultdict(lambda: defaultdict(lambda: defaultdict(float)))
    for condition, counts in sample_counts.items():
        for (locus_key, position), count in counts.items():
            positions_by_locus[locus_key][position][condition] += count

    joint_records = []
    excluded_missing_depth = 0
    excluded_not_apa = 0
    for locus_key, position_counts in positions_by_locus.items():
        clusters = cluster_positions(position_counts, max_pas_cluster_span)
        collapsed = []
        for members in clusters:
            pooled_weights = [sum(position_counts[position].values()) for position in members]
            representative = weighted_median(members, pooled_weights)
            counts = {
                condition: sum(position_counts[position].get(condition, 0.0) for position in members)
                for condition in sample_counts
            }
            collapsed.append((representative, members, counts))
        totals = {
            condition: sum(record[2][condition] for record in collapsed)
            for condition in sample_counts
        }
        if any(total < min_condition_locus_count for total in totals.values()):
            excluded_missing_depth += 1
            continue
        supported = 0
        for _representative, _members, counts in collapsed:
            if any(counts[condition] / totals[condition] >= min_isoform_fraction for condition in totals):
                supported += 1
        if supported < 2:
            excluded_not_apa += 1
            continue
        for representative, members, counts in collapsed:
            joint_records.append((locus_key, representative, members, counts, totals))

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
            ["isoform_id", "gene", "chrom", "strand", "anchor", "representative_pas", "member_pas_positions", "pooled_count"]
        )
        for locus_key, pas, members, counts, totals in sorted(joint_records):
            locus = lookup[locus_key]
            strand_label = "plus" if locus.strand == "+" else "minus"
            isoform_id = f"{locus.gene}__{locus.chrom}_{strand_label}_{locus.anchor}__PAS{pas}"
            start0, end = (locus.anchor, pas) if locus.strand == "+" else (pas, locus.anchor)
            pooled_count = sum(counts.values())
            bed_writer.writerow(
                [locus.chrom, start0, end, isoform_id, f"{pooled_count:g}", locus.strand, locus.gene]
            )
            for condition in sample_counts:
                usage_writer.writerow(
                    [isoform_id, condition, f"{counts[condition] / totals[condition]:.12g}"]
                )
            map_writer.writerow(
                [
                    isoform_id,
                    locus.gene,
                    locus.chrom,
                    locus.strand,
                    locus.anchor,
                    pas,
                    ",".join(map(str, members)),
                    f"{pooled_count:g}",
                ]
            )

    retained_loci = len({record[0] for record in joint_records})
    report = {
        "conditions": list(sample_counts),
        "sample_mapping_qc": sample_qc,
        "max_pas_cluster_span": max_pas_cluster_span,
        "min_cluster_count": min_cluster_count,
        "min_condition_locus_count": min_condition_locus_count,
        "min_isoform_fraction_for_apa": min_isoform_fraction,
        "excluded_loci_missing_condition_depth": excluded_missing_depth,
        "excluded_loci_without_two_supported_pas": excluded_not_apa,
        "retained_joint_loci": retained_loci,
        "retained_joint_isoforms": len(joint_records),
    }
    qc_out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", action="append", required=True, help="CONDITION=BED_PATH")
    parser.add_argument("--refgene", required=True, type=Path)
    parser.add_argument("--bed-out", required=True, type=Path)
    parser.add_argument("--usage-out", required=True, type=Path)
    parser.add_argument("--map-out", required=True, type=Path)
    parser.add_argument("--qc-out", required=True, type=Path)
    parser.add_argument("--max-pas-cluster-span", type=int, default=30)
    args = parser.parse_args()
    samples = []
    for value in args.sample:
        condition, separator, path = value.partition("=")
        if not separator or not condition or not path:
            parser.error(f"invalid --sample value: {value!r}")
        samples.append((condition, Path(path)))
    report = build_joint(
        samples,
        args.refgene,
        args.bed_out,
        args.usage_out,
        args.map_out,
        args.qc_out,
        max_pas_cluster_span=args.max_pas_cluster_span,
    )
    print(json.dumps(report, sort_keys=True))


if __name__ == "__main__":
    main()
