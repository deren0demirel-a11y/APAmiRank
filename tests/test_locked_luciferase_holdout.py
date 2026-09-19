from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

import pandas as pd


SCRIPT = Path(__file__).parents[1] / "benchmarks" / "locked_luciferase_holdout" / "lock_candidate_manifest.py"
SPEC = spec_from_file_location("lock_candidate_manifest", SCRIPT)
MODULE = module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_selection_is_deterministic_rank_independent_and_excludes_prior_pairs():
    rows = []
    for gene in ["A", "B", "C", "D"]:
        rows.append({
            "miRNA": "hsa-miR-16-5p", "gene": gene, "mti_ids": f"MIRT-{gene}",
            "reporter_assay": 1, "western_blot": 1, "qpcr": 0,
            "tier_a_or_b": True, "max_papers": 1,
            "primary_reporter_plus_wb_or_qpcr": True,
            "forbidden_rank_column": 1.0,
        })
    evidence = pd.DataFrame(rows)
    prior = pd.DataFrame({"miRNA": ["hsa-miR-16-5p"], "gene": ["A"]})
    first = MODULE.select_candidates(evidence, prior, per_mirna=2)
    second = MODULE.select_candidates(evidence.sample(frac=1, random_state=9), prior, per_mirna=2)
    assert first[["miRNA", "gene", "selection_hash"]].equals(
        second[["miRNA", "gene", "selection_hash"]]
    )
    assert "A" not in set(first.gene)
    assert "forbidden_rank_column" not in first.columns


def test_ineligible_evidence_is_not_selected():
    evidence = pd.DataFrame([{
        "miRNA": "hsa-miR-16-5p", "gene": "A", "mti_ids": "MIRT-A",
        "reporter_assay": 1, "western_blot": 0, "qpcr": 0,
        "tier_a_or_b": True, "max_papers": 1,
        "primary_reporter_plus_wb_or_qpcr": False,
    }])
    prior = pd.DataFrame(columns=["miRNA", "gene"])
    assert MODULE.select_candidates(evidence, prior).empty
