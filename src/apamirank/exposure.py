from __future__ import annotations

import itertools
from collections import defaultdict

from .io import read_tsv, write_tsv


REQUIRED_USAGE = {"isoform_id", "condition", "usage_fraction"}
REQUIRED_METADATA = {"isoform_id", "gene", "locus_id"}


def _fraction(value, field="usage_fraction"):
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Invalid numeric value in {field}: {value!r}") from exc
    if not 0.0 <= number <= 1.0:
        raise ValueError(f"{field} must be between 0 and 1; observed {number}")
    return number


def _load_usage(metadata_path, usage_path, required_loci, sum_tolerance=0.02):
    metadata = read_tsv(metadata_path)
    usage = read_tsv(usage_path)
    if not metadata:
        raise ValueError("Isoform metadata is empty")
    if not usage:
        raise ValueError("Isoform usage table is empty")
    missing_meta = REQUIRED_METADATA - set(metadata[0])
    missing_usage = REQUIRED_USAGE - set(usage[0])
    if missing_meta:
        raise ValueError("Missing metadata columns: " + ", ".join(sorted(missing_meta)))
    if missing_usage:
        raise ValueError("Missing usage columns: " + ", ".join(sorted(missing_usage)))

    isoform_info = {}
    locus_isoforms = defaultdict(set)
    for row in metadata:
        isoform = row["isoform_id"]
        if isoform in isoform_info:
            raise ValueError(f"Duplicate isoform_id in metadata: {isoform}")
        key = (row["gene"], row["locus_id"])
        isoform_info[isoform] = key
        locus_isoforms[key].add(isoform)

    by_context_locus = defaultdict(dict)
    conditions = set()
    for row in usage:
        isoform = row["isoform_id"]
        condition = row["condition"].strip()
        if not condition:
            raise ValueError("Usage condition cannot be empty")
        if isoform not in isoform_info:
            raise ValueError(f"Usage isoform missing from metadata: {isoform}")
        locus = isoform_info[isoform]
        key = (condition, locus)
        if isoform in by_context_locus[key]:
            raise ValueError(f"Duplicate usage row for {condition}/{isoform}")
        by_context_locus[key][isoform] = _fraction(row["usage_fraction"])
        conditions.add(condition)

    for condition in sorted(conditions):
        for locus in sorted(required_loci):
            expected = locus_isoforms.get(locus)
            if not expected:
                raise ValueError(f"Site locus absent from metadata: {locus[0]}/{locus[1]}")
            observed = set(by_context_locus.get((condition, locus), {}))
            missing = expected - observed
            if missing:
                raise ValueError(
                    f"Incomplete isoform usage for {condition}/{locus[0]}/{locus[1]}; "
                    f"missing {sorted(missing)}"
                )
            total = sum(by_context_locus[(condition, locus)].values())
            if abs(total - 1.0) > sum_tolerance:
                raise ValueError(
                    f"Usage fractions for {condition}/{locus[0]}/{locus[1]} sum to "
                    f"{total:.6f}, outside 1 +/- {sum_tolerance}"
                )
    return sorted(conditions), by_context_locus


def score_condition_exposure(
    site_ranking_path,
    metadata_path,
    usage_path,
    site_output_path,
    gene_output_path,
    contrast_output_path,
    sum_tolerance=0.02,
):
    """Weight ranked sites by observed isoform usage in each condition.

    The underlying binding-evidence score is descriptive, not a probability.
    Accordingly, the gene-level cumulative score is reported as an additive
    evidence summary and must not be interpreted as a probability of targeting.
    """
    sites = read_tsv(site_ranking_path)
    if not sites:
        raise ValueError("Site ranking is empty")
    required_site = {
        "site_id", "miRNA_id", "gene", "locus_id", "source_isoforms",
        "binding_evidence_score_equal_weight", "equal_weight_rank",
    }
    missing_site = required_site - set(sites[0])
    if missing_site:
        raise ValueError("Missing site-ranking columns: " + ", ".join(sorted(missing_site)))

    required_loci = {(row["gene"], row["locus_id"]) for row in sites}
    conditions, usage = _load_usage(
        metadata_path, usage_path, required_loci, sum_tolerance=float(sum_tolerance)
    )

    condition_sites = []
    for condition in conditions:
        rows = []
        for site in sites:
            locus = (site["gene"], site["locus_id"])
            carrying = {value for value in site["source_isoforms"].split(";") if value}
            exposure = sum(usage[(condition, locus)].get(isoform, 0.0) for isoform in carrying)
            evidence = float(site["binding_evidence_score_equal_weight"])
            row = dict(site)
            row.update({
                "condition": condition,
                "site_exposure_fraction": f"{exposure:.6f}",
                "exposure_weighted_evidence_score": f"{exposure * evidence:.6f}",
            })
            rows.append(row)
        rows.sort(
            key=lambda row: (
                -float(row["exposure_weighted_evidence_score"]),
                int(row["equal_weight_rank"]),
                row["site_id"],
            )
        )
        for rank, row in enumerate(rows, 1):
            row["condition_site_rank"] = rank
        condition_sites.extend(rows)
    write_tsv(site_output_path, condition_sites)

    by_condition_gene = defaultdict(list)
    for row in condition_sites:
        by_condition_gene[(row["condition"], row["gene"])].append(row)

    condition_genes = []
    for (condition, gene), rows in by_condition_gene.items():
        weighted = [float(row["exposure_weighted_evidence_score"]) for row in rows]
        exposures = [float(row["site_exposure_fraction"]) for row in rows]
        best = max(
            rows,
            key=lambda row: (
                float(row["exposure_weighted_evidence_score"]),
                -int(row["equal_weight_rank"]),
            ),
        )
        condition_genes.append({
            "condition": condition,
            "miRNA_id": best["miRNA_id"],
            "gene": gene,
            "number_canonical_sites": len(rows),
            "number_exposed_sites": sum(value > 0 for value in exposures),
            "best_site_id": best["site_id"],
            "best_site_exposure_fraction": best["site_exposure_fraction"],
            "best_exposure_weighted_evidence_score": best["exposure_weighted_evidence_score"],
            "cumulative_exposure_weighted_evidence_score": f"{sum(weighted):.6f}",
            "mean_site_exposure_fraction": f"{sum(exposures) / len(exposures):.6f}",
            "score_interpretation": "descriptive_additive_not_probability",
        })

    for condition in conditions:
        rows = [row for row in condition_genes if row["condition"] == condition]
        rows.sort(
            key=lambda row: (
                -float(row["cumulative_exposure_weighted_evidence_score"]),
                -float(row["best_exposure_weighted_evidence_score"]),
                row["gene"],
            )
        )
        for rank, row in enumerate(rows, 1):
            row["condition_gene_rank"] = rank
    condition_genes.sort(key=lambda row: (row["condition"], int(row["condition_gene_rank"])))
    write_tsv(gene_output_path, condition_genes)

    lookup = {
        (row["condition"], row["gene"]): row
        for row in condition_genes
    }
    genes = sorted({row["gene"] for row in condition_genes})
    contrasts = []
    for condition_a, condition_b in itertools.combinations(conditions, 2):
        for gene in genes:
            row_a = lookup.get((condition_a, gene))
            row_b = lookup.get((condition_b, gene))
            if not (row_a and row_b):
                continue
            score_a = float(row_a["cumulative_exposure_weighted_evidence_score"])
            score_b = float(row_b["cumulative_exposure_weighted_evidence_score"])
            contrasts.append({
                "miRNA_id": row_a["miRNA_id"],
                "gene": gene,
                "condition_a": condition_a,
                "condition_b": condition_b,
                "cumulative_score_a": f"{score_a:.6f}",
                "cumulative_score_b": f"{score_b:.6f}",
                "delta_cumulative_score_b_minus_a": f"{score_b - score_a:.6f}",
                "gene_rank_a": row_a["condition_gene_rank"],
                "gene_rank_b": row_b["condition_gene_rank"],
                "delta_rank_b_minus_a": int(row_b["condition_gene_rank"]) - int(row_a["condition_gene_rank"]),
            })
    contrasts.sort(
        key=lambda row: (
            row["condition_a"], row["condition_b"],
            -abs(float(row["delta_cumulative_score_b_minus_a"])), row["gene"],
        )
    )
    write_tsv(contrast_output_path, contrasts)
    return condition_sites, condition_genes, contrasts
