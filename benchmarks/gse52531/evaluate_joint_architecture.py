#!/usr/bin/env python3
"""Compare cognate and non-cognate AIR profiles on the common PAS atlas."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from evaluate_architecture import (
    CELL_FILE_LABELS,
    MIRNA_LABELS,
    analysis_role,
    load_accession_gene_map,
    load_expression,
    pearson,
    spearman,
)


def load_joint_scores(run_dir):
    unweighted = {}
    site_sums = {}
    with (run_dir / "02_site_ranking.tsv").open(encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            gene = row["gene"]
            site_sums[gene] = site_sums.get(gene, 0.0) + float(
                row["binding_evidence_score_equal_weight"]
            )
    unweighted.update(site_sums)
    weighted = {}
    with (run_dir / "03c_condition_gene_ranking.tsv").open(encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            weighted[(row["condition"], row["gene"])] = float(
                row["cumulative_exposure_weighted_evidence_score"]
            )
    return unweighted, weighted


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--refgene", required=True, type=Path)
    parser.add_argument("--expression-dir", required=True, type=Path)
    parser.add_argument("--runs-root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    accession_gene = load_accession_gene_map(args.refgene)
    results = []
    for short_mirna, mirna_label in MIRNA_LABELS.items():
        unweighted, weighted = load_joint_scores(args.runs_root / short_mirna)
        for outcome_cell in CELL_FILE_LABELS:
            expression_path = (
                args.expression_dir
                / f"GSE52530_{CELL_FILE_LABELS[outcome_cell]}.expData.qn.txt.gz"
            )
            expression = load_expression(expression_path, accession_gene, mirna_label)
            base_genes = {
                gene
                for gene in set(expression) & set(unweighted)
                if expression[gene]["control_mean_rpkm"] >= 1
                and expression[gene]["control_mean_rpkm"] > 0
                and expression[gene]["treatment_mean_rpkm"] > 0
            }
            observed = {
                gene: -__import__("math").log2(
                    expression[gene]["treatment_mean_rpkm"]
                    / expression[gene]["control_mean_rpkm"]
                )
                for gene in base_genes
            }
            for exposure_source in ["unweighted", *CELL_FILE_LABELS]:
                if exposure_source == "unweighted":
                    genes = sorted(base_genes)
                    scores = [unweighted[gene] for gene in genes]
                else:
                    genes = sorted(
                        gene for gene in base_genes if (exposure_source, gene) in weighted
                    )
                    scores = [weighted[(exposure_source, gene)] for gene in genes]
                values = [observed[gene] for gene in genes]
                r = pearson(scores, values)
                rho = spearman(scores, values)
                results.append(
                    {
                        "outcome_cell_line": outcome_cell,
                        "miRNA": mirna_label,
                        "analysis_role": analysis_role(outcome_cell, mirna_label),
                        "exposure_source": exposure_source,
                        "is_cognate_exposure": exposure_source == outcome_cell,
                        "n_genes": len(genes),
                        "pearson_r": r,
                        "pearson_r2": r**2,
                        "spearman_rho": rho,
                    }
                )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(results[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(results)
    print(args.output)


if __name__ == "__main__":
    main()
