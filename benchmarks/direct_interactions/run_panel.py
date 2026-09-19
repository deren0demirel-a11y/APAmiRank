#!/usr/bin/env python3
"""Run a predeclared five-miRNA direct-interaction benchmark panel."""

from __future__ import annotations

import argparse
from pathlib import Path

import yaml

from apamirank.cli import run_config


MIRNAS = {
    "hsa-miR-16-5p_MIMAT0000069": "UAGCAGCACGUAAAUAUUGGCG",
    "hsa-miR-17-5p_MIMAT0000070": "CAAAGUGCUUACAGUGCAGGUAG",
    "hsa-miR-25-3p_MIMAT0000081": "CAUUGCACUUGUCUCGGUCUGA",
    "hsa-miR-92a-3p_MIMAT0000092": "UAUUGCACUUGUCCCGGCCUGU",
    "hsa-miR-186-5p_MIMAT0000456": "CAAAGAAUUCUCCUUUUGGGCU",
}
CELL_LINES = ("HEK293", "Huh7")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--derived-dir", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    args = parser.parse_args()
    configs = args.output_root / "configs"
    configs.mkdir(parents=True, exist_ok=True)
    for cell in CELL_LINES:
        for mirna_id, sequence in MIRNAS.items():
            short_id = mirna_id.split("_", 1)[0]
            run_name = f"{cell}__{short_id}"
            config = {
                "miRNA": {"id": mirna_id, "sequence": sequence},
                "input": {
                    "isoform_metadata": str((args.derived_dir / f"{cell}.isoforms.tsv").resolve()),
                    "isoform_fasta": str((args.derived_dir / f"{cell}.isoforms.fa").resolve()),
                    "isoform_usage": str((args.derived_dir / f"{cell}.usage.tsv").resolve()),
                    "eligible_genes": None,
                    "external_site_annotations": None,
                    "posthoc_gene_annotations": None,
                    "thermodynamic_evidence": None,
                },
                "output_dir": str((args.output_root / run_name).resolve()),
                "flank_nt": 50,
                "usage_sum_tolerance": 0.02,
                "scoring": {
                    "mode": "architecture",
                    "bootstrap_iterations": 250,
                    "opportunity_iterations": 1000,
                    "random_seed": 20260908,
                },
            }
            config_path = configs / f"{run_name}.yaml"
            config_path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
            run_config(config_path)
            print(f"Completed {run_name}")


if __name__ == "__main__":
    main()
