from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

import pandas as pd
import pytest


SCRIPT = Path(__file__).parents[1] / "benchmarks" / "strict_luciferase" / "evaluate_strict_luciferase.py"
SPEC = spec_from_file_location("evaluate_strict_luciferase", SCRIPT)
MODULE = module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def audit_rows():
    return pd.DataFrame([
        {"miRNA": "hsa-miR-1", "gene": "A", "target_region": "3UTR", "wt_repressed": "yes", "mutant_rescue": "yes", "mirna_specificity": "exact", "verification_level": "full_text_verified", "eligibility": "strict", "primary_url": "https://example.org/a"},
        {"miRNA": "hsa-miR-2", "gene": "B", "target_region": "3UTR", "wt_repressed": "yes", "mutant_rescue": "yes", "mirna_specificity": "exact", "verification_level": "construct_corroborated", "eligibility": "corroborated_strict", "primary_url": "https://example.org/b"},
        {"miRNA": "hsa-miR-3", "gene": "C", "target_region": "3UTR", "wt_repressed": "yes", "mutant_rescue": "unclear", "mirna_specificity": "exact", "verification_level": "abstract_only", "eligibility": "excluded", "primary_url": "https://example.org/c"},
    ])


def test_scopes_are_separated_and_descriptive():
    ranks = pd.DataFrame({
        "miRNA": ["hsa-miR-1", "hsa-miR-2"],
        "gene": ["A", "B"],
        "consensus_apa_exposure_percentile": [0.95, 0.85],
        "consensus_unweighted_additive_percentile": [0.75, 0.90],
    })
    details, summary = MODULE.evaluate(audit_rows(), ranks)
    assert len(details[details.scope.eq("full_text_strict")]) == 1
    assert len(details[details.scope.eq("plus_construct_corroborated")]) == 2
    primary = summary[summary.scope.eq("full_text_strict")].iloc[0]
    assert primary.top_decile_pairs == 1
    assert primary.inferential_test == "not_run_post_selection_small_n"


def test_eligible_rows_must_have_mutant_rescue():
    audit = audit_rows()
    audit.loc[0, "mutant_rescue"] = "unclear"
    with pytest.raises(ValueError, match="mutant_rescue"):
        MODULE.validate_audit(audit)
