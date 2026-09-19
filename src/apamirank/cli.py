from __future__ import annotations

import argparse
import json
from pathlib import Path

import yaml

from . import __version__
from .aggregate import aggregate_genes
from .annotate import annotate_rankings
from .exposure import score_condition_exposure
from .io import sha256
from .prepare import prepare_qapa
from .scan import scan_sites
from .score import score_sites
from .thermo import run_thermodynamics


def _resolve(base, value):
    if value in {None, "", False}:
        return None
    path = Path(value)
    return path if path.is_absolute() else (base / path).resolve()


def run_config(config_path):
    config_path = Path(config_path).resolve()
    with open(config_path, encoding="utf-8") as handle:
        config = yaml.safe_load(handle)
    base = config_path.parent
    output_dir = _resolve(base, config.get("output_dir", "../results"))
    output_dir.mkdir(parents=True, exist_ok=True)
    inputs = config["input"]
    scoring = config.get("scoring", {})
    paths = {
        "sites": output_dir / "01_canonical_sites.tsv",
        "site_rank": output_dir / "02_site_ranking.tsv",
        "lofo": output_dir / "02_leave_one_feature_out.tsv",
        "gene_rank": output_dir / "03_gene_ranking.tsv",
        "site_final": output_dir / "04_site_ranking_annotated.tsv",
        "gene_final": output_dir / "04_gene_ranking_annotated.tsv",
    }
    usage_path = _resolve(base, inputs.get("isoform_usage"))
    if usage_path:
        paths.update({
            "condition_sites": output_dir / "03b_condition_site_ranking.tsv",
            "condition_genes": output_dir / "03c_condition_gene_ranking.tsv",
            "condition_contrasts": output_dir / "03d_condition_gene_contrasts.tsv",
        })
    scan_sites(
        _resolve(base, inputs["isoform_metadata"]), _resolve(base, inputs["isoform_fasta"]),
        config["miRNA"]["id"], config["miRNA"]["sequence"], paths["sites"],
        _resolve(base, inputs.get("eligible_genes")), flank=int(config.get("flank_nt", 50)),
    )
    mode = scoring.get("mode", "architecture")
    if mode == "thermodynamic":
        evidence = _resolve(base, inputs.get("thermodynamic_evidence"))
        if not evidence:
            evidence = output_dir / "01b_thermodynamic_evidence.tsv"
            tools = config.get("thermodynamics", {})
            run_thermodynamics(paths["sites"], config["miRNA"]["sequence"], evidence,
                               rnahybrid=tools.get("rnahybrid_executable", "RNAhybrid"),
                               rnaup=tools.get("rnaup_executable", "RNAup"),
                               threads=int(tools.get("threads", 1)), strict=bool(tools.get("strict", True)))
        scoring_input = evidence
    else:
        scoring_input = paths["sites"]
    score_sites(scoring_input, paths["site_rank"], paths["lofo"], mode=mode,
                bootstrap=int(scoring.get("bootstrap_iterations", 1000)),
                random_seed=int(scoring.get("random_seed", 1)))
    aggregate_genes(paths["site_rank"], paths["gene_rank"],
                    iterations=int(scoring.get("opportunity_iterations", 10000)),
                    random_seed=int(scoring.get("random_seed", 1)))
    if usage_path:
        score_condition_exposure(
            paths["site_rank"],
            _resolve(base, inputs["isoform_metadata"]),
            usage_path,
            paths["condition_sites"],
            paths["condition_genes"],
            paths["condition_contrasts"],
            sum_tolerance=float(config.get("usage_sum_tolerance", 0.02)),
        )
    annotate_rankings(paths["site_rank"], paths["gene_rank"], paths["site_final"], paths["gene_final"],
                      _resolve(base, inputs.get("external_site_annotations")),
                      _resolve(base, inputs.get("posthoc_gene_annotations")))
    manifest = {name: {"path": str(path), "sha256": sha256(path)} for name, path in paths.items()}
    with open(output_dir / "run_manifest.json", "w", encoding="utf-8") as handle:
        json.dump({"APAmiRank_version": __version__, "config": str(config_path), "outputs": manifest}, handle, indent=2)
    return paths


def build_parser():
    parser = argparse.ArgumentParser(prog="apamirank", description="APA-aware miRNA target ranking")
    parser.add_argument("--version", action="version", version=__version__)
    sub = parser.add_subparsers(dest="command", required=True)
    run = sub.add_parser("run", help="Run the configured workflow")
    run.add_argument("--config", required=True)
    prep = sub.add_parser("prepare-qapa", help="Extract strand-correct isoform sequences from QAPA BED and genome FASTA")
    prep.add_argument("--qapa-bed", required=True)
    prep.add_argument("--genome-fasta", required=True)
    prep.add_argument("--metadata-out", required=True)
    prep.add_argument("--fasta-out", required=True)
    scan = sub.add_parser("scan", help="Scan canonical sites")
    scan.add_argument("--metadata", required=True); scan.add_argument("--fasta", required=True)
    scan.add_argument("--mirna-id", required=True); scan.add_argument("--mirna-sequence", required=True)
    scan.add_argument("--output", required=True); scan.add_argument("--eligible-genes"); scan.add_argument("--flank", type=int, default=50)
    thermo = sub.add_parser("thermo", help="Run RNAhybrid and RNAup for scanned sites")
    thermo.add_argument("--sites", required=True); thermo.add_argument("--mirna-sequence", required=True)
    thermo.add_argument("--output", required=True); thermo.add_argument("--rnahybrid", default="RNAhybrid")
    thermo.add_argument("--rnaup", default="RNAup"); thermo.add_argument("--threads", type=int, default=1)
    thermo.add_argument("--allow-parse-failures", action="store_true")
    exposure = sub.add_parser("exposure", help="Weight ranked sites by condition-specific isoform usage")
    exposure.add_argument("--site-ranking", required=True)
    exposure.add_argument("--metadata", required=True)
    exposure.add_argument("--usage", required=True)
    exposure.add_argument("--site-output", required=True)
    exposure.add_argument("--gene-output", required=True)
    exposure.add_argument("--contrast-output", required=True)
    exposure.add_argument("--sum-tolerance", type=float, default=0.02)
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    if args.command == "run":
        run_config(args.config)
    elif args.command == "prepare-qapa":
        prepare_qapa(args.qapa_bed, args.genome_fasta, args.metadata_out, args.fasta_out)
    elif args.command == "scan":
        scan_sites(args.metadata, args.fasta, args.mirna_id, args.mirna_sequence, args.output, args.eligible_genes, args.flank)
    elif args.command == "thermo":
        run_thermodynamics(args.sites, args.mirna_sequence, args.output, args.rnahybrid, args.rnaup, args.threads, not args.allow_parse_failures)
    elif args.command == "exposure":
        score_condition_exposure(
            args.site_ranking, args.metadata, args.usage,
            args.site_output, args.gene_output, args.contrast_output,
            args.sum_tolerance,
        )


if __name__ == "__main__":
    main()
