from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

import pandas as pd
import pytest


SCRIPT = Path(__file__).parents[1] / "benchmarks" / "locked_luciferase_holdout" / "evaluate_locked_holdout.py"
SPEC = spec_from_file_location("evaluate_locked_holdout", SCRIPT)
MODULE = module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_pending_screening_blocks_rank_join(tmp_path):
    keys = {"miRNA": ["hsa-miR-1"], "gene": ["A"], "selection_hash": ["abc"]}
    manifest = pd.DataFrame(keys)
    screening = pd.DataFrame({
        **keys, "screening_status": ["pending"],
        "eligibility_decision": ["pending"],
    })
    path = tmp_path / "screening.tsv"
    screening.to_csv(path, sep="\t", index=False)
    with pytest.raises(ValueError, match="incomplete"):
        MODULE.validate_lock(manifest, screening, path, MODULE.sha256(path))


def test_checksum_mismatch_blocks_evaluation(tmp_path):
    keys = {"miRNA": ["hsa-miR-1"], "gene": ["A"], "selection_hash": ["abc"]}
    manifest = pd.DataFrame(keys)
    screening = pd.DataFrame({
        **keys, "screening_status": ["complete"],
        "eligibility_decision": ["exclude"], "exclusion_code": ["NO_MUTANT"],
    })
    path = tmp_path / "screening.tsv"
    screening.to_csv(path, sep="\t", index=False)
    with pytest.raises(ValueError, match="SHA-256"):
        MODULE.validate_lock(manifest, screening, path, "0" * 64)


def test_complete_conforming_screening_opens_gate(tmp_path):
    keys = {"miRNA": ["hsa-miR-1"], "gene": ["A"], "selection_hash": ["abc"]}
    manifest = pd.DataFrame(keys)
    screening = pd.DataFrame({
        **keys,
        "screening_status": ["complete"],
        "eligibility_decision": ["include"],
        "article_title": ["Original study"],
        "primary_url": ["https://example.org/original"],
        "species": ["Homo sapiens"],
        "mature_mirna_exact": ["yes"],
        "target_region": ["3UTR"],
        "reporter_vector": ["luciferase"],
        "wt_repression": ["yes"],
        "mutant_design": ["cognate-site mutation"],
        "mutant_rescue": ["yes"],
        "evidence_note": ["WT repression was lost in the mutant."],
        "reviewer": ["reviewer"],
        "reviewed_on": ["2026-09-10"],
        "exclusion_code": [""],
    })
    path = tmp_path / "screening.tsv"
    screening.to_csv(path, sep="\t", index=False)
    MODULE.validate_lock(manifest, screening, path, MODULE.sha256(path))


def test_included_record_requires_canonical_biological_fields(tmp_path):
    keys = {"miRNA": ["hsa-miR-1"], "gene": ["A"], "selection_hash": ["abc"]}
    manifest = pd.DataFrame(keys)
    screening = pd.DataFrame({
        **keys,
        "screening_status": ["complete"],
        "eligibility_decision": ["include"],
        "article_title": ["Original study"],
        "primary_url": ["https://example.org/original"],
        "species": ["Mus musculus"],
        "mature_mirna_exact": ["yes"],
        "target_region": ["3UTR"],
        "reporter_vector": ["luciferase"],
        "wt_repression": ["yes"],
        "mutant_design": ["cognate-site mutation"],
        "mutant_rescue": ["yes"],
        "evidence_note": ["WT repression was lost in the mutant."],
        "reviewer": ["reviewer"],
        "reviewed_on": ["2026-09-10"],
        "exclusion_code": [""],
    })
    path = tmp_path / "screening.tsv"
    screening.to_csv(path, sep="\t", index=False)
    with pytest.raises(ValueError, match="species=Homo sapiens"):
        MODULE.validate_lock(manifest, screening, path, MODULE.sha256(path))
