from pathlib import Path

import pytest

from apamirank.exposure import score_condition_exposure
from apamirank.io import read_tsv
from apamirank.scan import scan_sites
from apamirank.score import score_sites


ROOT = Path(__file__).resolve().parents[1]


def _ranked_sites(tmp_path):
    sites = tmp_path / "sites.tsv"
    ranking = tmp_path / "ranking.tsv"
    scan_sites(
        ROOT / "examples/minimal/isoforms.tsv",
        ROOT / "examples/minimal/isoforms.fa",
        "hsa-miR-example",
        "AUGCAUGCAUGCAUGCAUGCAU",
        sites,
    )
    score_sites(sites, ranking, tmp_path / "lofo.tsv", bootstrap=20, random_seed=3)
    return ranking


def test_condition_exposure_tracks_distal_isoform_usage(tmp_path):
    ranking = _ranked_sites(tmp_path)
    site_out = tmp_path / "condition_sites.tsv"
    gene_out = tmp_path / "condition_genes.tsv"
    contrast_out = tmp_path / "contrasts.tsv"
    score_condition_exposure(
        ranking,
        ROOT / "examples/minimal/isoforms.tsv",
        ROOT / "examples/minimal/isoform_usage.tsv",
        site_out,
        gene_out,
        contrast_out,
    )

    rows = read_tsv(site_out)
    distal_site = next(row for row in rows if row["gene"] == "GENEA" and row["extension_specific_yes_no"] == "Yes" and row["condition"] == "distal_context")
    proximal_site = next(row for row in rows if row["gene"] == "GENEA" and row["extension_specific_yes_no"] == "Yes" and row["condition"] == "proximal_context")
    assert float(distal_site["site_exposure_fraction"]) == pytest.approx(0.85)
    assert float(proximal_site["site_exposure_fraction"]) == pytest.approx(0.10)

    common = [row for row in rows if row["extension_specific_yes_no"] == "No"]
    assert common
    assert all(float(row["site_exposure_fraction"]) == pytest.approx(1.0) for row in common)

    genes = read_tsv(gene_out)
    gene_a_distal = next(row for row in genes if row["gene"] == "GENEA" and row["condition"] == "distal_context")
    gene_a_proximal = next(row for row in genes if row["gene"] == "GENEA" and row["condition"] == "proximal_context")
    assert float(gene_a_distal["cumulative_exposure_weighted_evidence_score"]) > float(gene_a_proximal["cumulative_exposure_weighted_evidence_score"])
    assert read_tsv(contrast_out)


def test_incomplete_usage_is_rejected(tmp_path):
    ranking = _ranked_sites(tmp_path)
    usage = tmp_path / "incomplete.tsv"
    usage.write_text(
        "isoform_id\tcondition\tusage_fraction\n"
        "GENEA_I1\tcondition_a\t1.0\n"
        "GENEB_I1\tcondition_a\t0.5\n"
        "GENEB_I2\tcondition_a\t0.5\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="Incomplete isoform usage"):
        score_condition_exposure(
            ranking,
            ROOT / "examples/minimal/isoforms.tsv",
            usage,
            tmp_path / "site.tsv",
            tmp_path / "gene.tsv",
            tmp_path / "contrast.tsv",
        )
