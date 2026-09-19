from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

import numpy as np
import pandas as pd


SCRIPT = Path(__file__).parents[1] / "benchmarks" / "locked_luciferase_holdout" / "evaluate_holdout_uncertainty.py"
SPEC = spec_from_file_location("evaluate_holdout_uncertainty", SCRIPT)
MODULE = module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_rank_uniform_null_is_deterministic_and_bounded():
    first = MODULE.simulate_rank_uniform_null({"mir-a": 2, "mir-b": 1}, {"mir-a": 11, "mir-b": 7}, 20, 4)
    second = MODULE.simulate_rank_uniform_null({"mir-a": 2, "mir-b": 1}, {"mir-a": 11, "mir-b": 7}, 20, 4)
    assert first.equals(second)
    assert first[MODULE.METRICS].ge(0).all().all()
    assert first[MODULE.METRICS].le(1).all().all()


def test_cluster_bootstrap_resamples_complete_cluster_summaries():
    clusters = pd.DataFrame({
        "miRNA": ["mir-a", "mir-b"],
        "observable_pairs": [2, 1],
        "mean_apa": [0.8, 0.6],
        "median_apa": [0.8, 0.6],
        "top_decile_fraction": [0.5, 0.0],
        "top_quintile_fraction": [1.0, 0.0],
        "median_apa_minus_unweighted": [0.1, -0.1],
    })
    draws = MODULE.bootstrap_clusters(clusters, iterations=25, seed=3)
    assert len(draws) == 25
    assert set(draws.mean_apa.unique()).issubset({0.6, 0.7, 0.8})


def test_holm_adjustment_is_monotone_in_sorted_p_values():
    raw = pd.Series([0.01, 0.03, np.nan, 0.02], index=list("abcd"))
    adjusted = MODULE.holm_adjust(raw)
    ordered = adjusted.loc[raw.dropna().sort_values().index].to_numpy()
    assert np.all(np.diff(ordered) >= 0)
    assert adjusted.isna().sum() == 1
