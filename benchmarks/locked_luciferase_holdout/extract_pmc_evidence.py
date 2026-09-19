#!/usr/bin/env python3
"""Extract reporter/mutant evidence passages from PMC full-text XML.

This is a screening aid only.  It deliberately does not assign eligibility.
"""

from __future__ import annotations

import argparse
import csv
import re
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path


EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
REPORTER_RE = re.compile(r"\b(luciferase|reporter|psiCHECK|pMIR|pmirGLO)\b", re.I)
DESIGN_RE = re.compile(
    r"\b(mutant|mutated|mutation|mutagen|delet(?:e|ed|ion)|wild[ -]?type|3[′'’-]?UTR|untranslated region)\b",
    re.I,
)


def normalized_text(element: ET.Element) -> str:
    return re.sub(r"\s+", " ", "".join(element.itertext())).strip()


def evidence_units(article: ET.Element) -> list[tuple[str, str]]:
    """Return de-duplicated relevant paragraphs/captions with section labels."""
    parent = {child: node for node in article.iter() for child in node}
    seen: set[str] = set()
    units: list[tuple[str, str]] = []
    for node in article.iter():
        tag = node.tag.rsplit("}", 1)[-1]
        if tag not in {"p", "caption"}:
            continue
        text = normalized_text(node)
        if len(text) < 40 or not REPORTER_RE.search(text):
            continue
        if not DESIGN_RE.search(text):
            continue
        if text in seen:
            continue
        seen.add(text)
        label = ""
        cur = parent.get(node)
        while cur is not None:
            if cur.tag.rsplit("}", 1)[-1] == "sec":
                title = cur.find("title")
                if title is not None:
                    label = normalized_text(title)
                    break
            cur = parent.get(cur)
        units.append((label, text))
    return units


def fetch_xml(pmcids: list[str], *, timeout: int = 60) -> bytes:
    params = urllib.parse.urlencode(
        {"db": "pmc", "id": ",".join(pmcids), "rettype": "full", "retmode": "xml"}
    )
    req = urllib.request.Request(
        f"{EFETCH}?{params}",
        headers={"User-Agent": "APAmiRank-literature-screening/0.2.18"},
    )
    with urllib.request.urlopen(req, timeout=timeout) as response:
        return response.read()


def article_id(article: ET.Element, kind: str) -> str:
    for node in article.findall(".//article-id"):
        if node.attrib.get("pub-id-type") == kind:
            return normalized_text(node)
    return ""


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidates", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--chunk-size", type=int, default=20)
    args = parser.parse_args()

    with args.candidates.open(newline="", encoding="utf-8") as handle:
        candidates = list(csv.DictReader(handle, delimiter="\t"))
    by_pmcid: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in candidates:
        if row.get("pmcid"):
            by_pmcid[row["pmcid"]].append(row)

    extracted: dict[str, tuple[str, str, list[tuple[str, str]]]] = {}
    ids = sorted(by_pmcid)
    for start in range(0, len(ids), args.chunk_size):
        chunk = ids[start : start + args.chunk_size]
        root = ET.fromstring(fetch_xml(chunk))
        for article in root.findall(".//article"):
            pmcid = article_id(article, "pmc")
            if pmcid and not pmcid.startswith("PMC"):
                pmcid = f"PMC{pmcid}"
            extracted[pmcid] = (
                article_id(article, "pmid"),
                article_id(article, "doi"),
                evidence_units(article),
            )
        if start + args.chunk_size < len(ids):
            time.sleep(0.4)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "miRNA", "gene", "mti_ids", "pmid", "pmcid", "doi",
        "article_title", "section", "evidence_excerpt",
    ]
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for pmcid in ids:
            pmid, doi, units = extracted.get(pmcid, ("", "", []))
            for candidate in by_pmcid[pmcid]:
                for section, excerpt in units or [("", "")]:
                    writer.writerow(
                        {
                            "miRNA": candidate["miRNA"],
                            "gene": candidate["gene"],
                            "mti_ids": candidate["mti_ids"],
                            "pmid": pmid or candidate["pmid"],
                            "pmcid": pmcid,
                            "doi": doi or candidate["doi"],
                            "article_title": candidate["article_title"],
                            "section": section,
                            "evidence_excerpt": excerpt,
                        }
                    )

    print(f"{len(ids)} PMC articles inspected; {sum(len(v[2]) for v in extracted.values())} evidence passages")


if __name__ == "__main__":
    main()
