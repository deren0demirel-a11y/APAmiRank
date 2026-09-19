#!/usr/bin/env python3
"""Run both miRNAs on the common GSE52527 PAS atlas."""

from __future__ import annotations

import argparse
from pathlib import Path

import yaml

from apamirank.cli import run_config
from run_architecture_benchmark import MIRNAS


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--joint-dir", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    args = parser.parse_args()
    args.output_root.mkdir(parents=True, exist_ok=True)
    for mirna_id, sequence in MIRNAS.items():
        short_id = mirna_id.split("_", 1)[0]
        output_dir = (args.output_root / short_id).resolve()
        config = {
            "miRNA": {"id": mirna_id, "sequence": sequence},
            "input": {
                "isoform_metadata": str((args.joint_dir / "joint.isoforms.tsv").resolve()),
                "isoform_fasta": str((args.joint_dir / "joint.isoforms.fa").resolve()),
                "isoform_usage": str((args.joint_dir / "joint.usage.tsv").resolve()),
                "eligible_genes": None,
                "external_site_annotations": None,
                "posthoc_gene_annotations": None,
                "thermodynamic_evidence": None,
            },
            "output_dir": str(output_dir),
            "flank_nt": 50,
            "usage_sum_tolerance": 0.02,
            "scoring": {
                "mode": "architecture",
                "bootstrap_iterations": 250,
                "opportunity_iterations": 1000,
                "random_seed": 20260907,
            },
        }
        config_path = args.output_root / f"{short_id}.yaml"
        config_path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
        run_config(config_path)
        print(f"Completed joint {short_id}")


if __name__ == "__main__":
    main()
