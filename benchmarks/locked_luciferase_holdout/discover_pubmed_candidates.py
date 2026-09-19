#!/usr/bin/env python3
"""Discover PubMed candidates for locked pairs without making eligibility calls."""

from __future__ import annotations

import argparse
import json
import re
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

import pandas as pd


EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
ALIASES = {
    "hsa-miR-124-3p": ["miR-124", "miR-124-3p", "microRNA-124"],
    "hsa-miR-155-5p": ["miR-155", "miR-155-5p", "microRNA-155"],
    "hsa-miR-16-5p": ["miR-16", "miR-16-5p", "microRNA-16"],
    "hsa-miR-17-5p": ["miR-17", "miR-17-5p", "microRNA-17"],
    "hsa-miR-186-5p": ["miR-186", "miR-186-5p", "microRNA-186"],
    "hsa-miR-25-3p": ["miR-25", "miR-25-3p", "microRNA-25"],
    "hsa-miR-92a-3p": ["miR-92a", "miR-92a-3p", "microRNA-92a"],
}


def get(url: str) -> bytes:
    request = urllib.request.Request(url, headers={"User-Agent": "APAmiRank-literature-audit/0.2.18"})
    with urllib.request.urlopen(request, timeout=45) as response:
        return response.read()


def post(endpoint: str, parameters: dict) -> bytes:
    payload = urllib.parse.urlencode(parameters).encode("utf-8")
    request = urllib.request.Request(
        f"{EUTILS}/{endpoint}", data=payload,
        headers={"User-Agent": "APAmiRank-literature-audit/0.2.18"},
    )
    with urllib.request.urlopen(request, timeout=60) as response:
        return response.read()


def build_query(mirna: str, gene: str) -> str:
    aliases = ALIASES[mirna]
    mirna_clause = " OR ".join(f'"{alias}"[Title/Abstract]' for alias in aliases)
    return f"({mirna_clause}) AND \"{gene}\"[Title/Abstract]"


def search(query: str, retmax: int = 1000) -> list[str]:
    payload = json.loads(post("esearch.fcgi", {
        "db": "pubmed", "term": query, "retmax": retmax, "retmode": "json",
    }))
    return payload["esearchresult"]["idlist"]


def text_content(element: ET.Element | None) -> str:
    return "" if element is None else "".join(element.itertext()).strip()


def parse_articles(xml_bytes: bytes) -> dict[str, dict]:
    root = ET.fromstring(xml_bytes)
    records = {}
    for article in root.findall(".//PubmedArticle"):
        citation = article.find("MedlineCitation")
        pmid = text_content(citation.find("PMID"))
        record = citation.find("Article")
        title = text_content(record.find("ArticleTitle"))
        abstract = " ".join(text_content(x) for x in record.findall("Abstract/AbstractText"))
        journal = text_content(record.find("Journal/Title"))
        year = text_content(record.find("Journal/JournalIssue/PubDate/Year"))
        if not year:
            medline_date = text_content(record.find("Journal/JournalIssue/PubDate/MedlineDate"))
            match = re.search(r"(?:19|20)\d{2}", medline_date)
            year = match.group(0) if match else ""
        ids = {
            node.attrib.get("IdType", ""): text_content(node)
            for node in article.findall("PubmedData/ArticleIdList/ArticleId")
        }
        publication_types = ";".join(
            text_content(node) for node in record.findall("PublicationTypeList/PublicationType")
        )
        records[pmid] = {
            "pmid": pmid, "pmcid": ids.get("pmc", ""), "doi": ids.get("doi", ""),
            "article_title": title, "abstract": abstract, "journal": journal,
            "publication_year": year, "publication_types": publication_types,
        }
    return records


def priority_score(record: dict, mirna: str, gene: str) -> int:
    title = record["article_title"].lower()
    abstract = record["abstract"].lower()
    aliases = [x.lower() for x in ALIASES[mirna]]
    score = 0
    if gene.lower() in title:
        score += 3
    if any(alias in title for alias in aliases):
        score += 3
    if "luciferase" in abstract:
        score += 4
    if any(term in abstract for term in ["3′utr", "3'utr", "3′-utr", "3'-utr", "3' untranslated"]):
        score += 2
    if any(term in abstract for term in ["mutant", "mutation", "mutated"]):
        score += 2
    if "target" in abstract:
        score += 1
    return score


def record_matches_pair(record: dict, mirna: str, gene: str) -> bool:
    text = f"{record['article_title']} {record['abstract']}".lower()
    gene_match = re.search(rf"(?<![a-z0-9]){re.escape(gene.lower())}(?![a-z0-9])", text)
    mirna_match = any(alias.lower() in text for alias in ALIASES[mirna])
    return bool(gene_match and mirna_match)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--query-log", required=True, type=Path)
    args = parser.parse_args()
    manifest = pd.read_csv(args.manifest, sep="\t", dtype=str)
    clauses = [build_query(row.miRNA, row.gene) for row in manifest.itertuples(index=False)]
    union_query = " OR ".join(f"({clause})" for clause in clauses)
    all_pmids = set(search(union_query))

    records = {}
    ordered_pmids = sorted(all_pmids, key=int)
    for start in range(0, len(ordered_pmids), 100):
        ids = ",".join(ordered_pmids[start:start + 100])
        records.update(parse_articles(post("efetch.fcgi", {
            "db": "pubmed", "id": ids, "retmode": "xml",
        })))

    output_rows, query_rows = [], []
    for row, clause in zip(manifest.itertuples(index=False), clauses):
        matches = [record for record in records.values() if record_matches_pair(record, row.miRNA, row.gene)]
        matches.sort(key=lambda record: (-priority_score(record, row.miRNA, row.gene), -int(record["pmid"])))
        query_rows.append({
            "miRNA": row.miRNA, "gene": row.gene, "pair_query_clause": clause,
            "union_pubmed_hits": len(all_pmids), "locally_mapped_pair_hits": len(matches),
        })
        for position, record in enumerate(matches, start=1):
            output_rows.append({
                "miRNA": row.miRNA, "gene": row.gene, "mti_ids": row.mti_ids,
                "pubmed_query_position": position,
                "discovery_priority_score": priority_score(record, row.miRNA, row.gene),
                **record,
                "pubmed_url": f"https://pubmed.ncbi.nlm.nih.gov/{record['pmid']}/",
            })
    output = pd.DataFrame(output_rows)
    if not output.empty:
        output = output.sort_values(
            ["miRNA", "gene", "discovery_priority_score", "pubmed_query_position"],
            ascending=[True, True, False, True],
        )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(args.output, sep="\t", index=False)
    pd.DataFrame(query_rows).to_csv(args.query_log, sep="\t", index=False)
    print(f"{len(manifest)} locked pairs; {len(output)} pair-article candidates; {len(all_pmids)} unique PMIDs")


if __name__ == "__main__":
    main()
