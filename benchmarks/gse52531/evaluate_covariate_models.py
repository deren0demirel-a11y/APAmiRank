#!/usr/bin/env python3
"""Evaluate APAmiRank scores after adjustment for prespecified gene covariates.

The primary comparison is outcome-blind and uses the same genes for all models:

    C       = baseline expression + UTR length + site count + mean local AU
    C + U   = C + unweighted APAmiRank architecture score
    C + W   = C + exposure-weighted APAmiRank score
    C + U+W = C + both APAmiRank scores

Repeated K-fold predictions are pooled within each repeat before R-squared is
calculated.  The reported interval is the empirical 2.5--97.5 percentile range
across repeats; it is a stability interval, not a confidence interval.
"""

from __future__ import annotations

import argparse
import csv
import math
import random
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy import stats


BASE_COVARIATES = (
    "log10_control_rpkm",
    "log1p_max_benchmark_utr_length",
    "log1p_canonical_site_count",
    "mean_site_au_fraction",
)
MODELS = {
    "covariates_only": BASE_COVARIATES,
    "covariates_plus_unweighted": BASE_COVARIATES + ("unweighted_score",),
    "covariates_plus_weighted": BASE_COVARIATES + ("weighted_score",),
    "covariates_plus_both": BASE_COVARIATES
    + ("unweighted_score", "weighted_score"),
}


def percentile(values, probability):
    return float(np.quantile(np.asarray(values, dtype=float), probability))


def r_squared(observed, predicted):
    observed = np.asarray(observed, dtype=float)
    predicted = np.asarray(predicted, dtype=float)
    total = float(np.sum((observed - observed.mean()) ** 2))
    return 1.0 - float(np.sum((observed - predicted) ** 2)) / total if total > 0 else math.nan


def load_gene_architecture(run_dir: Path):
    summary = defaultdict(lambda: {"site_count": 0, "au_sum": 0.0})
    with (run_dir / "02_site_ranking.tsv").open(encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            record = summary[row["gene"]]
            record["site_count"] += 1
            record["au_sum"] += float(row["AU_flank30"])
    return {
        gene: {
            "canonical_site_count": values["site_count"],
            "mean_site_au_fraction": values["au_sum"] / values["site_count"],
        }
        for gene, values in summary.items()
    }


def load_max_utr_lengths(metadata: Path):
    lengths = defaultdict(int)
    with metadata.open(encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            length = int(row["end"]) - int(row["start0"])
            lengths[row["gene"]] = max(lengths[row["gene"]], length)
    return dict(lengths)


def load_evaluation_rows(path: Path, architecture, utr_lengths):
    rows = []
    with path.open(encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            gene = row["gene"]
            if gene not in architecture or gene not in utr_lengths:
                continue
            control = float(row["control_mean_rpkm"])
            if control <= 0:
                continue
            record = dict(row)
            record.update(architecture[gene])
            record["max_benchmark_utr_length"] = utr_lengths[gene]
            record["log10_control_rpkm"] = math.log10(control)
            record["log1p_max_benchmark_utr_length"] = math.log1p(utr_lengths[gene])
            record["log1p_canonical_site_count"] = math.log1p(
                architecture[gene]["canonical_site_count"]
            )
            for key in ("repression", "unweighted_score", "weighted_score"):
                record[key] = float(record[key])
            rows.append(record)
    return rows


def standardized_design(rows, features):
    matrix = np.asarray([[float(row[key]) for key in features] for row in rows], dtype=float)
    center = matrix.mean(axis=0)
    scale = matrix.std(axis=0)
    scale[scale == 0] = 1.0
    return (matrix - center) / scale


def fit_ols(rows, features):
    y = np.asarray([row["repression"] for row in rows], dtype=float)
    x = standardized_design(rows, features)
    design = np.column_stack([np.ones(len(rows)), x])
    beta, _, _, _ = np.linalg.lstsq(design, y, rcond=None)
    fitted = design @ beta
    residual = y - fitted
    rss = float(residual @ residual)
    r2 = r_squared(y, fitted)
    p = design.shape[1] - 1
    adjusted = 1 - (1 - r2) * (len(rows) - 1) / (len(rows) - p - 1)
    covariance = np.linalg.pinv(design.T @ design) * (rss / (len(rows) - p - 1))
    standard_errors = np.sqrt(np.diag(covariance))
    t_values = np.divide(beta, standard_errors, out=np.zeros_like(beta), where=standard_errors > 0)
    p_values = 2 * stats.t.sf(np.abs(t_values), df=len(rows) - p - 1)
    return {
        "rss": rss,
        "r2": r2,
        "adjusted_r2": adjusted,
        "beta": beta,
        "p_values": p_values,
        "condition_number": float(np.linalg.cond(design)),
    }


def nested_f_test(reduced, full, n, p_reduced, p_full):
    numerator_df = p_full - p_reduced
    denominator_df = n - p_full - 1
    improvement = max(0.0, reduced["rss"] - full["rss"])
    statistic = (improvement / numerator_df) / (full["rss"] / denominator_df)
    return statistic, float(stats.f.sf(statistic, numerator_df, denominator_df))


def repeated_kfold(rows, repeats=100, folds=10, seed=20260908):
    y = np.asarray([row["repression"] for row in rows], dtype=float)
    raw = {
        name: np.asarray([[float(row[key]) for key in features] for row in rows], dtype=float)
        for name, features in MODELS.items()
    }
    by_model = defaultdict(list)
    deltas = defaultdict(list)
    for repeat in range(repeats):
        indices = list(range(len(rows)))
        random.Random(seed + repeat).shuffle(indices)
        fold_ids = np.empty(len(rows), dtype=int)
        for position, index in enumerate(indices):
            fold_ids[index] = position % folds
        predictions = {name: np.empty(len(rows), dtype=float) for name in MODELS}
        for fold in range(folds):
            test = fold_ids == fold
            train = ~test
            for name, matrix in raw.items():
                center = matrix[train].mean(axis=0)
                scale = matrix[train].std(axis=0)
                scale[scale == 0] = 1.0
                x_train = (matrix[train] - center) / scale
                x_test = (matrix[test] - center) / scale
                design_train = np.column_stack([np.ones(train.sum()), x_train])
                design_test = np.column_stack([np.ones(test.sum()), x_test])
                beta, _, _, _ = np.linalg.lstsq(design_train, y[train], rcond=None)
                predictions[name][test] = design_test @ beta
        scores = {name: r_squared(y, prediction) for name, prediction in predictions.items()}
        for name, score in scores.items():
            by_model[name].append(score)
        deltas["weighted_minus_unweighted"].append(
            scores["covariates_plus_weighted"] - scores["covariates_plus_unweighted"]
        )
        deltas["weighted_beyond_unweighted"].append(
            scores["covariates_plus_both"] - scores["covariates_plus_unweighted"]
        )
        deltas["unweighted_beyond_weighted"].append(
            scores["covariates_plus_both"] - scores["covariates_plus_weighted"]
        )
    return by_model, deltas


def write_gene_table(path: Path, rows):
    fields = [
        "cell_line", "miRNA", "analysis_role", "gene", "selected_refseq_id",
        "control_mean_rpkm", "treatment_mean_rpkm", "repression",
        "unweighted_score", "weighted_score", "canonical_site_count",
        "mean_site_au_fraction", "max_benchmark_utr_length",
        "log10_control_rpkm", "log1p_max_benchmark_utr_length",
        "log1p_canonical_site_count",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows({key: row[key] for key in fields} for row in rows)


def evaluate(evaluation_dir: Path, runs_root: Path, derived_dir: Path, output_dir: Path,
             repeats=100, folds=10, seed=20260908):
    output_dir.mkdir(parents=True, exist_ok=True)
    summary = []
    for evaluation_path in sorted(evaluation_dir.glob("*.gene_evaluation.tsv")):
        run_name = evaluation_path.name.removesuffix(".gene_evaluation.tsv")
        cell_line = run_name.split("__", 1)[0]
        architecture = load_gene_architecture(runs_root / run_name)
        utr_lengths = load_max_utr_lengths(derived_dir / f"{cell_line}.isoforms.tsv")
        all_rows = load_evaluation_rows(evaluation_path, architecture, utr_lengths)
        for threshold in (0.1, 1.0, 5.0, 10.0):
            rows = [row for row in all_rows if float(row["control_mean_rpkm"]) >= threshold]
            if len(rows) <= len(MODELS["covariates_plus_both"]) + 2:
                continue
            if threshold == 1.0:
                write_gene_table(output_dir / f"{run_name}.covariate_input.tsv", rows)
            fits = {name: fit_ols(rows, features) for name, features in MODELS.items()}
            cv, deltas = repeated_kfold(rows, repeats=repeats, folds=folds, seed=seed)
            f_w, p_w = nested_f_test(
                fits["covariates_plus_unweighted"], fits["covariates_plus_both"],
                len(rows), len(MODELS["covariates_plus_unweighted"]),
                len(MODELS["covariates_plus_both"]),
            )
            weighted_index = 1 + MODELS["covariates_plus_both"].index("weighted_score")
            row = {
                "cell_line": rows[0]["cell_line"],
                "miRNA": rows[0]["miRNA"],
                "analysis_role": rows[0]["analysis_role"],
                "control_rpkm_threshold": threshold,
                "n_genes": len(rows),
                "cv_repeats": repeats,
                "cv_folds": folds,
            }
            for name in MODELS:
                row[f"in_sample_r2__{name}"] = fits[name]["r2"]
                row[f"adjusted_r2__{name}"] = fits[name]["adjusted_r2"]
                row[f"cv_r2_median__{name}"] = percentile(cv[name], 0.5)
                row[f"cv_r2_p025__{name}"] = percentile(cv[name], 0.025)
                row[f"cv_r2_p975__{name}"] = percentile(cv[name], 0.975)
            for name, values in deltas.items():
                row[f"cv_delta_median__{name}"] = percentile(values, 0.5)
                row[f"cv_delta_p025__{name}"] = percentile(values, 0.025)
                row[f"cv_delta_p975__{name}"] = percentile(values, 0.975)
                row[f"cv_delta_positive_fraction__{name}"] = sum(value > 0 for value in values) / len(values)
            row.update({
                "partial_f_weighted_beyond_unweighted": f_w,
                "partial_f_p_weighted_beyond_unweighted": p_w,
                "standardized_beta_weighted_in_full": fits["covariates_plus_both"]["beta"][weighted_index],
                "standardized_beta_weighted_p_in_full": fits["covariates_plus_both"]["p_values"][weighted_index],
                "full_model_condition_number": fits["covariates_plus_both"]["condition_number"],
            })
            summary.append(row)
    output = output_dir / "covariate_model_summary.tsv"
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(summary)
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evaluation-dir", required=True, type=Path)
    parser.add_argument("--runs-root", required=True, type=Path)
    parser.add_argument("--derived-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--cv-repeats", type=int, default=100)
    parser.add_argument("--cv-folds", type=int, default=10)
    parser.add_argument("--seed", type=int, default=20260908)
    args = parser.parse_args()
    print(evaluate(args.evaluation_dir, args.runs_root, args.derived_dir, args.output_dir,
                   repeats=args.cv_repeats, folds=args.cv_folds, seed=args.seed))


if __name__ == "__main__":
    main()
