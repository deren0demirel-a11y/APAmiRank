#!/usr/bin/env python3
"""Evaluate weighted versus unweighted APAmiRank scores against GSE52530."""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import random
from collections import defaultdict
from pathlib import Path


MIRNA_LABELS = {"hsa-miR-124-3p": "miR-124", "hsa-miR-155-5p": "miR-155"}
CELL_FILE_LABELS = {"HeLa": "HeLa", "HEK293": "HEK293", "Huh7": "HUH7", "IMR90": "IMR90"}


def mean(values):
    return sum(values) / len(values)


def pearson(x, y):
    if len(x) < 3:
        return float("nan")
    xbar, ybar = mean(x), mean(y)
    numerator = sum((a - xbar) * (b - ybar) for a, b in zip(x, y, strict=True))
    dx = sum((a - xbar) ** 2 for a in x)
    dy = sum((b - ybar) ** 2 for b in y)
    return numerator / math.sqrt(dx * dy) if dx > 0 and dy > 0 else float("nan")


def ranks(values):
    order = sorted(range(len(values)), key=lambda index: values[index])
    result = [0.0] * len(values)
    start = 0
    while start < len(order):
        end = start + 1
        while end < len(order) and values[order[end]] == values[order[start]]:
            end += 1
        average_rank = (start + 1 + end) / 2
        for position in range(start, end):
            result[order[position]] = average_rank
        start = end
    return result


def spearman(x, y):
    return pearson(ranks(x), ranks(y))


def percentile(values, probability):
    ordered = sorted(values)
    position = (len(ordered) - 1) * probability
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[lower]
    return ordered[lower] * (upper - position) + ordered[upper] * (position - lower)


def bootstrap_delta_r2(rows, iterations=2000, seed=20260907):
    generator = random.Random(seed)
    deltas = []
    for _ in range(iterations):
        sample = [rows[generator.randrange(len(rows))] for _ in rows]
        observed = [row["repression"] for row in sample]
        unweighted = [row["unweighted_score"] for row in sample]
        weighted = [row["weighted_score"] for row in sample]
        r_unweighted = pearson(unweighted, observed)
        r_weighted = pearson(weighted, observed)
        if math.isfinite(r_unweighted) and math.isfinite(r_weighted):
            deltas.append(r_weighted**2 - r_unweighted**2)
    return percentile(deltas, 0.025), percentile(deltas, 0.975)


def load_accession_gene_map(refgene: Path):
    genes = defaultdict(set)
    with gzip.open(refgene, "rt", encoding="utf-8") as handle:
        for row in csv.reader(handle, delimiter="\t"):
            genes[row[1]].add(row[12])
    return {accession: next(iter(values)) for accession, values in genes.items() if len(values) == 1}


def load_expression(expression_path: Path, accession_gene, mirna_label: str):
    transcripts = []
    with gzip.open(expression_path, "rt", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        groups = {column.split(",")[0].rsplit("_", 1)[0]: index for index, column in enumerate(header)}
        treatment_key = mirna_label
        control_key = "Puc19"
        if "Puc19_124" in groups:
            control_key = "Puc19_124" if mirna_label == "miR-124" else "Puc19_155"
        for row in reader:
            accession = row[0]
            gene = accession_gene.get(accession)
            if gene is None:
                continue
            control = [float(value) for value in row[groups[control_key]].split(",")]
            treatment = [float(value) for value in row[groups[treatment_key]].split(",")]
            transcripts.append(
                {
                    "gene": gene,
                    "refseq_id": accession,
                    "control_mean_rpkm": mean(control),
                    "treatment_mean_rpkm": mean(treatment),
                }
            )
    # Select the most highly expressed RefSeq record using control data only.
    selected = {}
    for row in transcripts:
        key = (row["control_mean_rpkm"], row["refseq_id"])
        old = selected.get(row["gene"])
        if old is None or key > (old["control_mean_rpkm"], old["refseq_id"]):
            selected[row["gene"]] = row
    return selected


def load_scores(run_dir: Path):
    unweighted = defaultdict(float)
    with (run_dir / "02_site_ranking.tsv").open(encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            unweighted[row["gene"]] += float(row["binding_evidence_score_equal_weight"])
    weighted = {}
    with (run_dir / "03c_condition_gene_ranking.tsv").open(encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            weighted[row["gene"]] = float(row["cumulative_exposure_weighted_evidence_score"])
    return unweighted, weighted


def analysis_role(cell_line, mirna_label):
    if cell_line in {"HeLa", "HEK293", "Huh7"}:
        return "primary"
    if cell_line == "IMR90" and mirna_label == "miR-124":
        return "supporting"
    return "exploratory_not_primary_paper_design"


def evaluate(refgene: Path, expression_dir: Path, runs_root: Path, output_dir: Path):
    accession_gene = load_accession_gene_map(refgene)
    output_dir.mkdir(parents=True, exist_ok=True)
    summary = []
    for cell_line in CELL_FILE_LABELS:
        expression_path = expression_dir / f"GSE52530_{CELL_FILE_LABELS[cell_line]}.expData.qn.txt.gz"
        for short_mirna, mirna_label in MIRNA_LABELS.items():
            run_dir = runs_root / f"{cell_line}__{short_mirna}"
            expression = load_expression(expression_path, accession_gene, mirna_label)
            unweighted, weighted = load_scores(run_dir)
            genes = sorted(set(expression) & set(unweighted) & set(weighted))
            rows = []
            for gene in genes:
                record = expression[gene]
                control = record["control_mean_rpkm"]
                treatment = record["treatment_mean_rpkm"]
                if control <= 0 or treatment <= 0:
                    continue
                log2fc = math.log2(treatment / control)
                rows.append(
                    {
                        "cell_line": cell_line,
                        "miRNA": mirna_label,
                        "analysis_role": analysis_role(cell_line, mirna_label),
                        "gene": gene,
                        "selected_refseq_id": record["refseq_id"],
                        "control_mean_rpkm": control,
                        "treatment_mean_rpkm": treatment,
                        "log2fc": log2fc,
                        "repression": -log2fc,
                        "unweighted_score": unweighted[gene],
                        "weighted_score": weighted[gene],
                    }
                )
            gene_path = output_dir / f"{cell_line}__{short_mirna}.gene_evaluation.tsv"
            with gene_path.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
                writer.writeheader()
                writer.writerows(rows)

            for threshold in (0.1, 1.0, 5.0, 10.0):
                eligible = [row for row in rows if row["control_mean_rpkm"] >= threshold]
                observed = [row["repression"] for row in eligible]
                score_u = [row["unweighted_score"] for row in eligible]
                score_w = [row["weighted_score"] for row in eligible]
                pearson_u, pearson_w = pearson(score_u, observed), pearson(score_w, observed)
                spearman_u, spearman_w = spearman(score_u, observed), spearman(score_w, observed)
                ci_low = ci_high = float("nan")
                if threshold == 1.0 and len(eligible) >= 3:
                    ci_low, ci_high = bootstrap_delta_r2(eligible)
                summary.append(
                    {
                        "cell_line": cell_line,
                        "miRNA": mirna_label,
                        "analysis_role": analysis_role(cell_line, mirna_label),
                        "control_rpkm_threshold": threshold,
                        "n_site_containing_genes": len(eligible),
                        "pearson_r_unweighted": pearson_u,
                        "pearson_r_weighted": pearson_w,
                        "pearson_r2_unweighted": pearson_u**2,
                        "pearson_r2_weighted": pearson_w**2,
                        "delta_pearson_r2_weighted_minus_unweighted": pearson_w**2 - pearson_u**2,
                        "spearman_rho_unweighted": spearman_u,
                        "spearman_rho_weighted": spearman_w,
                        "delta_r2_bootstrap_ci_low": ci_low,
                        "delta_r2_bootstrap_ci_high": ci_high,
                    }
                )
    summary_path = output_dir / "architecture_benchmark_summary.tsv"
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(summary)
    return summary_path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--refgene", required=True, type=Path)
    parser.add_argument("--expression-dir", required=True, type=Path)
    parser.add_argument("--runs-root", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    print(evaluate(args.refgene, args.expression_dir, args.runs_root, args.output_dir))


if __name__ == "__main__":
    main()
