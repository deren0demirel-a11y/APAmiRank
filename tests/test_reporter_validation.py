from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

import numpy as np
import pandas as pd


SCRIPT = Path(__file__).parents[1] / "benchmarks" / "reporter_validation" / "evaluate_reporter_validation.py"
SPEC = spec_from_file_location("evaluate_reporter_validation", SCRIPT)
MODULE = module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_auc_from_arrays_handles_separation_and_ties():
    assert MODULE.auc_from_arrays(np.array([0.8, 0.9]), np.array([0.1, 0.2])) == 1.0
    assert MODULE.auc_from_arrays(np.array([0.5]), np.array([0.5])) == 0.5


def test_small_positive_set_suppresses_bootstrap():
    ranks = pd.DataFrame({
        "gene": ["A", "B", "C", "D", "E", "F"],
        "consensus_apa_exposure_percentile": [1, .8, .6, .4, .2, 0],
        "consensus_unweighted_additive_percentile": [1, .8, .6, .4, .2, 0],
    })
    result = MODULE.bootstrap_auc_increment(ranks, {"A", "B"}, iterations=10, seed=1)
    assert np.isnan(result["apa_minus_unweighted_auc"])
