#!/usr/bin/env python3
"""Run APAmiRank architecture mode for the GSE52531 cell-line/miRNA grid."""

from __future__ import annotations

import argparse
from pathlib import Path

import yaml

from apamirank.cli import run_config


MIRNAS = {
    "hsa-miR-124-3p_MIMAT0000422": "UAAGGCACGCGGUGAAUGCCAA",
    "hsa-miR-155-5p_MIMAT0000646": "UUAAUGCUAAUCGUGAUAGGGGUU",
}
CELL_LINES = ("HeLa", "HEK293", "Huh7", "IMR90")


def run_grid(derived_dir: Path, output_root: Path):
    configs_dir = output_root / "configs"
    configs_dir.mkdir(parents=True, exist_ok=True)
    completed = []
    for cell_line in CELL_LINES:
        for mirna_id, sequence in MIRNAS.items():
            short_id = mirna_id.split("_", 1)[0]
            run_name = f"{cell_line}__{short_id}"
            config = {
                "miRNA": {"id": mirna_id, "sequence": sequence},
                "input": {
                    "isoform_metadata": str((derived_dir / f"{cell_line}.isoforms.tsv").resolve()),
                    "isoform_fasta": str((derived_dir / f"{cell_line}.isoforms.fa").resolve()),
                    "isoform_usage": str((derived_dir / f"{cell_line}.usage.tsv").resolve()),
                    "eligible_genes": None,
                    "external_site_annotations": None,
                    "posthoc_gene_annotations": None,
                    "thermodynamic_evidence": None,
                },
                "output_dir": str((output_root / run_name).resolve()),
                "flank_nt": 50,
                "usage_sum_tolerance": 0.02,
                "scoring": {
                    "mode": "architecture",
                    "bootstrap_iterations": 250,
                    "opportunity_iterations": 1000,
                    "random_seed": 20260907,
                },
            }
            config_path = configs_dir / f"{run_name}.yaml"
            config_path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
            run_config(config_path)
            completed.append(run_name)
            print(f"Completed {run_name}")
    return completed


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--derived-dir", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    args = parser.parse_args()
    run_grid(args.derived_dir, args.output_root)


if __name__ == "__main__":
    main()
