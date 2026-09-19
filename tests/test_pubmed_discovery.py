from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / "benchmarks" / "locked_luciferase_holdout" / "discover_pubmed_candidates.py"
SPEC = spec_from_file_location("discover_pubmed_candidates", SCRIPT)
MODULE = module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_query_is_fielded_and_pair_specific():
    query = MODULE.build_query("hsa-miR-124-3p", "GRB2")
    assert '"GRB2"[Title/Abstract]' in query
    assert '"miR-124"[Title/Abstract]' in query
    assert "luciferase" not in query.lower()


def test_priority_is_discovery_only_and_rewards_assay_terms():
    base = {"article_title": "miR-124 targets GRB2", "abstract": "Direct target."}
    assay = {"article_title": "miR-124 targets GRB2", "abstract": "Luciferase 3'-UTR mutant target assay."}
    assert MODULE.priority_score(assay, "hsa-miR-124-3p", "GRB2") > MODULE.priority_score(base, "hsa-miR-124-3p", "GRB2")


def test_pair_mapping_requires_both_gene_and_mirna():
    record = {"article_title": "miR-124 targets GRB2", "abstract": ""}
    assert MODULE.record_matches_pair(record, "hsa-miR-124-3p", "GRB2")
    assert not MODULE.record_matches_pair(record, "hsa-miR-124-3p", "MAPK14")
