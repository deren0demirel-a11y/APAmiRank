from __future__ import annotations

import re
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from .io import read_tsv, write_tsv
from .seed import normalize_mirna


MFE_RE = re.compile(r"^mfe:\s+(-?\d+(?:\.\d+)?)\s+kcal/mol", re.M)
PVAL_RE = re.compile(r"^p-value:\s+([0-9.eE+-]+)", re.M)
POS_RE = re.compile(r"^position\s+(\d+)\s*$", re.M)
UP_RE = re.compile(
    r"(?P<t1>\d+),(?P<t2>\d+)\s*:\s*(?P<q1>\d+),(?P<q2>\d+)\s*"
    r"\(\s*(?P<total>-?\d+(?:\.\d+)?)\s*=\s*"
    r"(?P<int>-?\d+(?:\.\d+)?)\s*\+\s*(?P<open>-?\d+(?:\.\d+)?)"
    r"(?:\s*\+\s*(?P<open_short>-?\d+(?:\.\d+)?))?\s*\)"
)


def _overlap(a1, a2, b1, b2):
    return max(0, min(a2, b2) - max(a1, b1) + 1)


def _nt_count(text):
    return sum(base in "ACGU" for base in text.upper())


def parse_rnahybrid(text, expected_start, expected_end):
    mfe, pvalue, position = MFE_RE.search(text), PVAL_RE.search(text), POS_RE.search(text)
    if not (mfe and position):
        return {"rnahybrid_parse_status": "failed"}
    lines = text.splitlines()
    position_line = next((i for i, line in enumerate(lines) if POS_RE.match(line.strip())), None)
    if position_line is None or position_line + 4 >= len(lines):
        return {"rnahybrid_parse_status": "failed_alignment"}
    span = _nt_count(lines[position_line + 1]) + _nt_count(lines[position_line + 2])
    start = int(position.group(1)); end = start + span - 1
    overlap = _overlap(start, end, expected_start, expected_end)
    expected_length = expected_end - expected_start + 1
    full = start <= expected_start and end >= expected_end
    concordance = "canonical_seed_fully_covered" if full else ("canonical_seed_partially_overlapped" if overlap else "canonical_seed_not_overlapped")
    return {
        "rnahybrid_parse_status": "ok", "rnahybrid_mfe_kcal_mol": mfe.group(1),
        "rnahybrid_p_value": pvalue.group(1) if pvalue else "",
        "rnahybrid_position_start1": start, "rnahybrid_position_end1": end,
        "rnahybrid_canonical_overlap_nt": overlap,
        "rnahybrid_canonical_overlap_fraction": f"{overlap / expected_length:.6f}",
        "rnahybrid_canonical_fully_covered": "Yes" if full else "No",
        "rnahybrid_coordinate_concordance_class": concordance,
    }


def parse_rnaup(text, expected_start, expected_end):
    matches = list(UP_RE.finditer(text))
    if not matches:
        return {"rnaup_parse_status": "failed"}
    values = matches[-1].groupdict()
    start, end = int(values["t1"]), int(values["t2"])
    overlap = _overlap(start, end, expected_start, expected_end)
    expected_length = expected_end - expected_start + 1
    return {
        "rnaup_parse_status": "ok", "rnaup_target_start1": start, "rnaup_target_end1": end,
        "rnaup_query_start1": values["q1"], "rnaup_query_end1": values["q2"],
        "rnaup_total_dG_kcal_mol": values["total"],
        "rnaup_interaction_dG_kcal_mol": values["int"],
        "rnaup_target_opening_dG_kcal_mol": values["open"],
        "rnaup_short_opening_dG_kcal_mol": values["open_short"] or "",
        "rnaup_canonical_overlap_nt": overlap,
        "rnaup_canonical_overlap_fraction": f"{overlap / expected_length:.6f}",
        "rnaup_canonical_fully_covered": "Yes" if start <= expected_start and end >= expected_end else "No",
        "rnaup_canonical_any_overlap": "Yes" if overlap else "No",
    }


def run_thermodynamics(sites_path, mirna_sequence, output_path, rnahybrid="RNAhybrid", rnaup="RNAup", threads=1, strict=True):
    sites = read_tsv(sites_path)
    mirna = normalize_mirna(mirna_sequence)
    if not sites:
        raise ValueError("Site table is empty")

    def evaluate(row):
        target = row["window_sequence_RNA_5to3"].upper().replace("T", "U")
        expected_start, expected_end = int(row["site_start1_window"]), int(row["site_end1_window"])
        with tempfile.TemporaryDirectory(prefix="apamirank_") as tmp:
            tmp = Path(tmp)
            target_fa, mirna_fa = tmp / "target.fa", tmp / "mirna.fa"
            target_fa.write_text(f">{row['site_id']}\n{target}\n", encoding="utf-8")
            mirna_fa.write_text(f">miRNA\n{mirna}\n", encoding="utf-8")
            hybrid = subprocess.run([rnahybrid, "-s", "3utr_human", "-b", "1", "-t", str(target_fa), "-q", str(mirna_fa)], text=True, capture_output=True)
            if hybrid.returncode:
                raise RuntimeError(f"RNAhybrid failed for {row['site_id']}: {hybrid.stderr.strip()}")
            up = subprocess.run([rnaup, "-o"], input=f">{row['site_id']}\n{target}&{mirna}\n", text=True, capture_output=True)
            if up.returncode:
                raise RuntimeError(f"RNAup failed for {row['site_id']}: {up.stderr.strip()}")
        result = dict(row)
        result.update(parse_rnahybrid(hybrid.stdout, expected_start, expected_end))
        result.update(parse_rnaup(up.stdout, expected_start, expected_end))
        return result

    with ThreadPoolExecutor(max_workers=max(1, int(threads))) as pool:
        output = list(pool.map(evaluate, sites))
    failures = [row["site_id"] for row in output if row.get("rnahybrid_parse_status") != "ok" or row.get("rnaup_parse_status") != "ok"]
    if strict and failures:
        raise RuntimeError(f"Thermodynamic output parsing failed for {len(failures)} sites; first: {failures[:5]}")
    write_tsv(output_path, output)
    return output
