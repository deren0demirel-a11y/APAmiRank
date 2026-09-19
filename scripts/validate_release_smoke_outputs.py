#!/usr/bin/env python3
"""Fail-closed validation of the APAmiRank thermodynamic smoke-test outputs."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path


EXPECTED_VERSION = "0.2.18"
REQUIRED_FILES = (
    "01_canonical_sites.tsv",
    "01b_thermodynamic_evidence.tsv",
    "02_site_ranking.tsv",
    "02_leave_one_feature_out.tsv",
    "03_gene_ranking.tsv",
    "03b_condition_site_ranking.tsv",
    "03c_condition_gene_ranking.tsv",
    "03d_condition_gene_contrasts.tsv",
    "04_site_ranking_annotated.tsv",
    "04_gene_ranking_annotated.tsv",
    "run_manifest.json",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def validate(results: Path) -> dict:
    missing = [name for name in REQUIRED_FILES if not (results / name).is_file()]
    if missing:
        raise RuntimeError(f"Missing required outputs: {missing}")

    manifest_path = results / "run_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("APAmiRank_version") != EXPECTED_VERSION:
        raise RuntimeError(
            f"Expected APAmiRank {EXPECTED_VERSION}, found {manifest.get('APAmiRank_version')!r}"
        )

    outputs = manifest.get("outputs", {})
    evidence_record = outputs.get("thermodynamic_evidence")
    if not evidence_record:
        raise RuntimeError("Manifest does not record thermodynamic_evidence")

    checksum_failures = []
    for name, record in outputs.items():
        path = Path(record["path"])
        if not path.is_file():
            checksum_failures.append(f"{name}: missing {path}")
            continue
        observed = sha256(path)
        if observed != record.get("sha256"):
            checksum_failures.append(f"{name}: expected {record.get('sha256')}, observed {observed}")
    if checksum_failures:
        raise RuntimeError("Manifest checksum failure(s): " + "; ".join(checksum_failures))

    evidence = read_tsv(results / "01b_thermodynamic_evidence.tsv")
    if not evidence:
        raise RuntimeError("Thermodynamic evidence table is empty")
    parse_failures = [
        row.get("site_id", "<unknown>")
        for row in evidence
        if row.get("rnahybrid_parse_status") != "ok" or row.get("rnaup_parse_status") != "ok"
    ]
    if parse_failures:
        raise RuntimeError(f"Thermodynamic parse failures: {parse_failures[:10]}")

    rankings = read_tsv(results / "02_site_ranking.tsv")
    genes = read_tsv(results / "03_gene_ranking.tsv")
    if len(rankings) != len(evidence):
        raise RuntimeError(
            f"Site-count mismatch: evidence={len(evidence)}, ranking={len(rankings)}"
        )
    if not genes:
        raise RuntimeError("Gene ranking is empty")

    return {
        "status": "PASS",
        "APAmiRank_version": EXPECTED_VERSION,
        "results_directory": str(results.resolve()),
        "site_count": len(rankings),
        "gene_count": len(genes),
        "rnahybrid_parse_failures": 0,
        "rnaup_parse_failures": 0,
        "manifest_outputs_verified": len(outputs),
        "manifest_sha256": sha256(manifest_path),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    args = parser.parse_args()
    report = validate(args.results.resolve())
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
