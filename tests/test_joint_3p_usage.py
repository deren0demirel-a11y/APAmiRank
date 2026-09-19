import importlib.util
import sys
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / "benchmarks" / "gse52531" / "prepare_joint_3p_usage.py"
sys.path.insert(0, str(SCRIPT.parent))
SPEC = importlib.util.spec_from_file_location("prepare_joint_3p_usage", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def test_pas_clustering_has_bounded_span_and_weighted_median():
    clusters = MODULE.cluster_positions({100: {}, 120: {}, 129: {}, 131: {}, 161: {}}, 30)
    assert clusters == [[100, 120, 129], [131, 161]]
    assert MODULE.weighted_median([100, 110, 120], [1, 8, 1]) == 110
