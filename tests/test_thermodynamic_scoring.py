import json
from pathlib import Path

import yaml

from apamirank.cli import run_config
from apamirank.io import read_tsv, write_tsv
from apamirank.scan import scan_sites
from apamirank.score import score_sites
from apamirank.thermo import parse_rnaup


ROOT = Path(__file__).resolve().parents[1]


def test_thermodynamic_scoring_schema_without_external_binaries(tmp_path):
    sites_path = tmp_path / "sites.tsv"
    rows = scan_sites(
        ROOT / "examples/minimal/isoforms.tsv",
        ROOT / "examples/minimal/isoforms.fa",
        "hsa-miR-example",
        "AUGCAUGCAUGCAUGCAUGCAU",
        sites_path,
    )
    for index, row in enumerate(rows):
        row.update({
            "rnahybrid_coordinate_concordance_class": "canonical_seed_fully_covered",
            "rnahybrid_mfe_kcal_mol": str(-18.0 - index),
            "rnaup_canonical_fully_covered": "Yes",
            "rnaup_canonical_any_overlap": "Yes",
            "rnaup_total_dG_kcal_mol": str(-12.0 - index),
            "rnaup_target_opening_dG_kcal_mol": str(3.0 + index),
        })
    evidence = tmp_path / "thermo.tsv"
    write_tsv(evidence, rows)
    ranked = tmp_path / "ranked.tsv"
    lofo = tmp_path / "lofo.tsv"
    score_sites(evidence, ranked, lofo, mode="thermodynamic", bootstrap=25, random_seed=4)
    output = read_tsv(ranked)
    assert all(row["binding_evidence_tier"] == "consensus_full" for row in output)
    assert all(len(row["features_in_primary_score"].split(";")) == 7 for row in output)


def test_rnaup_parser_extracts_energy_and_overlap():
    parsed = parse_rnaup("1,7 : 2,8 (-10.0 = -15.0 + 5.0)\n", 1, 7)
    assert parsed["rnaup_parse_status"] == "ok"
    assert parsed["rnaup_canonical_fully_covered"] == "Yes"
    assert parsed["rnaup_total_dG_kcal_mol"] == "-10.0"


def test_precomputed_thermodynamic_evidence_is_checksummed(tmp_path):
    sites_path = tmp_path / "sites.tsv"
    rows = scan_sites(
        ROOT / "examples/minimal/isoforms.tsv",
        ROOT / "examples/minimal/isoforms.fa",
        "hsa-miR-example",
        "AUGCAUGCAUGCAUGCAUGCAU",
        sites_path,
    )
    for index, row in enumerate(rows):
        row.update({
            "rnahybrid_parse_status": "ok",
            "rnahybrid_coordinate_concordance_class": "canonical_seed_fully_covered",
            "rnahybrid_mfe_kcal_mol": str(-18.0 - index),
            "rnaup_parse_status": "ok",
            "rnaup_canonical_fully_covered": "Yes",
            "rnaup_canonical_any_overlap": "Yes",
            "rnaup_total_dG_kcal_mol": str(-12.0 - index),
            "rnaup_target_opening_dG_kcal_mol": str(3.0 + index),
        })
    evidence = tmp_path / "thermodynamic_evidence.tsv"
    write_tsv(evidence, rows)
    config = {
        "miRNA": {"id": "hsa-miR-example", "sequence": "AUGCAUGCAUGCAUGCAUGCAU"},
        "input": {
            "isoform_metadata": str(ROOT / "examples/minimal/isoforms.tsv"),
            "isoform_fasta": str(ROOT / "examples/minimal/isoforms.fa"),
            "thermodynamic_evidence": str(evidence),
        },
        "output_dir": str(tmp_path / "results"),
        "scoring": {
            "mode": "thermodynamic",
            "bootstrap_iterations": 25,
            "opportunity_iterations": 50,
            "random_seed": 19,
        },
    }
    config_path = tmp_path / "config.yaml"
    config_path.write_text(yaml.safe_dump(config), encoding="utf-8")

    paths = run_config(config_path)
    manifest = json.loads((tmp_path / "results/run_manifest.json").read_text(encoding="utf-8"))

    assert paths["thermodynamic_evidence"] == evidence
    assert manifest["outputs"]["thermodynamic_evidence"]["path"] == str(evidence)
    assert len(manifest["outputs"]["thermodynamic_evidence"]["sha256"]) == 64
