import importlib.util
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / "benchmarks" / "gse52531" / "evaluate_covariate_models.py"
SPEC = importlib.util.spec_from_file_location("evaluate_covariate_models", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_repeated_cv_detects_added_weighted_signal():
    rows = []
    for index in range(120):
        covariate = (index % 17) / 17
        unweighted = (index % 11) / 11
        exposure = ((index * 7) % 19) / 19
        weighted = 0.55 * unweighted + exposure
        rows.append(
            {
                "repression": 0.2 * covariate + 0.5 * unweighted + 1.1 * exposure,
                "log10_control_rpkm": covariate,
                "log1p_max_benchmark_utr_length": (index % 13) / 13,
                "log1p_canonical_site_count": (index % 5) / 5,
                "mean_site_au_fraction": (index % 23) / 23,
                "unweighted_score": unweighted,
                "weighted_score": weighted,
            }
        )
    _, deltas = MODULE.repeated_kfold(rows, repeats=10, folds=5, seed=7)
    assert min(deltas["weighted_beyond_unweighted"]) > 0


def test_nested_f_test_returns_probability():
    reduced = {"rss": 20.0}
    full = {"rss": 10.0}
    statistic, probability = MODULE.nested_f_test(reduced, full, n=100, p_reduced=5, p_full=6)
    assert statistic > 0
    assert 0 <= probability <= 1
