from __future__ import annotations

import math
import random
from collections import defaultdict

from .io import read_tsv, write_tsv


def _quantile(values, probability):
    position = (len(values) - 1) * probability
    lower, upper = math.floor(position), math.ceil(position)
    if lower == upper:
        return values[lower]
    return values[lower] * (upper - position) + values[upper] * (position - lower)


def aggregate_genes(site_ranking, output_path, iterations=10000, random_seed=1):
    sites = read_tsv(site_ranking)
    if not sites:
        raise ValueError("Site ranking is empty")
    by_gene = defaultdict(list)
    all_scores = []
    for row in sites:
        score = float(row["binding_evidence_score_equal_weight"])
        all_scores.append(score)
        by_gene[row["gene"]].append((score, row))
    rng = random.Random(random_seed)
    null_by_count = {}
    for count in sorted({len(values) for values in by_gene.values()}):
        null_by_count[count] = sorted(max(rng.sample(all_scores, count)) for _ in range(iterations))
    output = []
    for gene, values in by_gene.items():
        best_score, best = max(values, key=lambda item: (item[0], -int(item[1]["equal_weight_rank"])))
        null = null_by_count[len(values)]
        tail = (sum(value >= best_score for value in null) + 1) / (iterations + 1)
        output.append({
            "miRNA_id": best["miRNA_id"],
            "gene": gene,
            "number_canonical_sites": len(values),
            "best_site_id": best["site_id"],
            "best_seed_class": best["seed_class"],
            "best_binding_evidence_tier": best["binding_evidence_tier"],
            "best_site_coordinate": best["genomic_coordinate_1based"],
            "best_equal_weight_score": f"{best_score:.6f}",
            "best_equal_weight_rank": best["equal_weight_rank"],
            "site_count_adjusted_empirical_percentile": f"{1 - tail:.6f}",
            "site_count_adjusted_empirical_tail_probability": f"{tail:.6g}",
            "matched_null_median_max_score": f"{_quantile(null, 0.5):.6f}",
            "matched_null_q95_max_score": f"{_quantile(null, 0.95):.6f}",
            "extension_specific_sites": sum(row["extension_specific_yes_no"] == "Yes" for _, row in values),
            "common_region_sites": sum(row["region_class"] == "common" for _, row in values),
            "score_model": best["score_model"],
        })
    output.sort(key=lambda row: (float(row["site_count_adjusted_empirical_tail_probability"]), -float(row["best_equal_weight_score"]), row["gene"]))
    for rank, row in enumerate(output, 1):
        row["opportunity_adjusted_descriptive_rank"] = rank
    write_tsv(output_path, output)
    return output
