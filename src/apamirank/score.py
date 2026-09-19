from __future__ import annotations

import math
import random
from array import array

from .io import read_tsv, write_tsv


THERMO_FEATURES = ["seed", "au", "hybrid_concordance", "hybrid_mfe", "rnaup_concordance", "rnaup_total", "rnaup_opening"]
ARCH_FEATURES = ["seed", "au"]


def _number(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def percentile(values, lower_better=False):
    good = sorted((value, index) for index, value in enumerate(values) if not math.isnan(value))
    output = [0.0] * len(values)
    n = len(good)
    if not n:
        return output
    cursor = 0
    while cursor < n:
        stop = cursor + 1
        while stop < n and good[stop][0] == good[cursor][0]:
            stop += 1
        p = ((cursor + stop - 1) / 2) / (n - 1) if n > 1 else 1.0
        if lower_better:
            p = 1 - p
        for j in range(cursor, stop):
            output[good[j][1]] = p
        cursor = stop
    return output


def ordinal_ranks(scores):
    order = sorted(range(len(scores)), key=lambda i: (-scores[i], i))
    result = [0] * len(scores)
    for rank, index in enumerate(order, 1):
        result[index] = rank
    return result


def _spearman(a, b):
    mean_a, mean_b = sum(a) / len(a), sum(b) / len(b)
    numerator = sum((x - mean_a) * (y - mean_b) for x, y in zip(a, b))
    denom_a = sum((x - mean_a) ** 2 for x in a)
    denom_b = sum((y - mean_b) ** 2 for y in b)
    return numerator / math.sqrt(denom_a * denom_b) if denom_a and denom_b else math.nan


def score_sites(input_path, output_path, lofo_path, mode="architecture", bootstrap=1000, random_seed=1):
    rows = read_tsv(input_path)
    if not rows:
        raise ValueError("Site table is empty; no canonical sites were detected")
    mode = mode.lower()
    if mode not in {"architecture", "thermodynamic"}:
        raise ValueError("Scoring mode must be architecture or thermodynamic")
    n = len(rows)
    seed = [(_number(row["seed_strength_ordinal"]) - 1) / 3 for row in rows]
    au = percentile([_number(row["AU_flank30"]) for row in rows])
    features = {"seed": seed, "au": au}
    feature_names = ARCH_FEATURES

    if mode == "thermodynamic":
        required = {"rnahybrid_coordinate_concordance_class", "rnahybrid_mfe_kcal_mol", "rnaup_canonical_fully_covered", "rnaup_canonical_any_overlap", "rnaup_total_dG_kcal_mol", "rnaup_target_opening_dG_kcal_mol"}
        missing = required - set(rows[0])
        if missing:
            raise ValueError("Thermodynamic evidence columns are missing: " + ", ".join(sorted(missing)))
        hybrid_concordance, hybrid_mfe, rnaup_concordance, rnaup_total, rnaup_opening = [], [], [], [], []
        for row in rows:
            h_class = row["rnahybrid_coordinate_concordance_class"]
            h_score = 1.0 if h_class == "canonical_seed_fully_covered" else (0.5 if h_class == "canonical_seed_partially_overlapped" else 0.0)
            hybrid_concordance.append(h_score)
            hybrid_mfe.append(_number(row["rnahybrid_mfe_kcal_mol"]) if h_score else math.nan)
            full = row["rnaup_canonical_fully_covered"] == "Yes"
            any_overlap = row["rnaup_canonical_any_overlap"] == "Yes"
            u_score = 1.0 if full else (0.5 if any_overlap else 0.0)
            rnaup_concordance.append(u_score)
            rnaup_total.append(_number(row["rnaup_total_dG_kcal_mol"]) if any_overlap else math.nan)
            rnaup_opening.append(_number(row["rnaup_target_opening_dG_kcal_mol"]) if any_overlap else math.nan)
        features.update({
            "hybrid_concordance": hybrid_concordance,
            "hybrid_mfe": percentile(hybrid_mfe, lower_better=True),
            "rnaup_concordance": rnaup_concordance,
            "rnaup_total": percentile(rnaup_total, lower_better=True),
            "rnaup_opening": percentile(rnaup_opening, lower_better=True),
        })
        feature_names = THERMO_FEATURES

    scores = [sum(features[name][i] for name in feature_names) / len(feature_names) for i in range(n)]
    base_ranks = ordinal_ranks(scores)
    top_n = max(1, math.ceil(n * 0.10))
    base_top = {i for i, rank in enumerate(base_ranks) if rank <= top_n}
    lofo = []
    if len(feature_names) > 1:
        for omitted in feature_names:
            retained = [name for name in feature_names if name != omitted]
            alternative = [sum(features[name][i] for name in retained) / len(retained) for i in range(n)]
            alt_ranks = ordinal_ranks(alternative)
            alt_top = {i for i, rank in enumerate(alt_ranks) if rank <= top_n}
            lofo.append({
                "omitted_feature": omitted,
                "spearman_rank_correlation": f"{_spearman(base_ranks, alt_ranks):.6f}",
                "top_decile_retained_fraction": f"{len(base_top & alt_top) / top_n:.6f}",
            })

    rng = random.Random(random_seed)
    bootstrap_ranks = [array("I") for _ in rows]
    bootstrap_top = [0] * n
    for _ in range(bootstrap):
        weights = [rng.expovariate(1) for _ in feature_names]
        total = sum(weights)
        weights = [weight / total for weight in weights]
        alternative = [sum(weights[j] * features[name][i] for j, name in enumerate(feature_names)) for i in range(n)]
        alt_ranks = ordinal_ranks(alternative)
        for i, rank in enumerate(alt_ranks):
            bootstrap_ranks[i].append(rank)
            bootstrap_top[i] += rank <= top_n

    output = []
    for i, row in enumerate(rows):
        sorted_ranks = sorted(bootstrap_ranks[i])
        entry = dict(row)
        entry.update({
            "score_model": mode,
            "binding_evidence_tier": "architecture_only",
            "features_in_primary_score": ";".join(feature_names),
            "binding_evidence_score_equal_weight": f"{scores[i]:.6f}",
            "equal_weight_rank": base_ranks[i],
            "bootstrap_rank_q05": sorted_ranks[int(0.05 * (bootstrap - 1))],
            "bootstrap_rank_median": sorted_ranks[int(0.50 * (bootstrap - 1))],
            "bootstrap_rank_q95": sorted_ranks[int(0.95 * (bootstrap - 1))],
            "bootstrap_top_decile_frequency": f"{bootstrap_top[i] / bootstrap:.6f}",
            "seed_score": f"{seed[i]:.6f}",
            "AU_context_percentile": f"{au[i]:.6f}",
        })
        if mode == "thermodynamic":
            hybrid = features["hybrid_concordance"][i]
            rnaup = features["rnaup_concordance"][i]
            if hybrid == 1 and rnaup == 1:
                entry["binding_evidence_tier"] = "consensus_full"
            elif hybrid > 0 and rnaup > 0:
                entry["binding_evidence_tier"] = "consensus_overlap"
            elif hybrid == 1 or rnaup == 1:
                entry["binding_evidence_tier"] = "single_tool_full"
            elif hybrid > 0 or rnaup > 0:
                entry["binding_evidence_tier"] = "single_tool_overlap"
            else:
                entry["binding_evidence_tier"] = "seed_only"
        for name in feature_names:
            if name not in {"seed", "au"}:
                entry[name + "_score"] = f"{features[name][i]:.6f}"
        output.append(entry)
    output.sort(key=lambda row: int(row["equal_weight_rank"]))
    write_tsv(output_path, output)
    write_tsv(lofo_path, lofo, ["omitted_feature", "spearman_rank_correlation", "top_decile_retained_fraction"])
    return output
