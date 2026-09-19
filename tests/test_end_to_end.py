from pathlib import Path

import yaml

from apamirank.cli import run_config
from apamirank.io import read_tsv


ROOT = Path(__file__).resolve().parents[1]


def test_architecture_workflow_and_posthoc_rank_separation(tmp_path):
    config = {
        "miRNA": {"id": "hsa-miR-example", "sequence": "AUGCAUGCAUGCAUGCAUGCAU"},
        "input": {
            "isoform_metadata": str(ROOT / "examples/minimal/isoforms.tsv"),
            "isoform_fasta": str(ROOT / "examples/minimal/isoforms.fa"),
            "isoform_usage": str(ROOT / "examples/minimal/isoform_usage.tsv"),
            "posthoc_gene_annotations": str(ROOT / "examples/minimal/posthoc_genes.tsv"),
        },
        "output_dir": str(tmp_path / "results"),
        "scoring": {"mode": "architecture", "bootstrap_iterations": 50, "opportunity_iterations": 100, "random_seed": 17},
    }
    config_path = tmp_path / "config.yaml"
    config_path.write_text(yaml.safe_dump(config), encoding="utf-8")
    paths = run_config(config_path)
    original = read_tsv(paths["gene_rank"])
    annotated = read_tsv(paths["gene_final"])
    assert [row["opportunity_adjusted_descriptive_rank"] for row in original] == [row["opportunity_adjusted_descriptive_rank"] for row in annotated]
    gene_a = next(row for row in annotated if row["gene"] == "GENEA")
    assert gene_a["posthoc_gene_annotation_present"] == "Yes"
    assert gene_a["posthoc_log2FoldChange"] == "-1.25"
    assert (tmp_path / "results/run_manifest.json").exists()
    assert (tmp_path / "results/03b_condition_site_ranking.tsv").exists()
    assert (tmp_path / "results/03c_condition_gene_ranking.tsv").exists()
    assert (tmp_path / "results/03d_condition_gene_contrasts.tsv").exists()
