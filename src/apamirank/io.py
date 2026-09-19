from __future__ import annotations

import csv
import gzip
import hashlib
from pathlib import Path
from typing import Iterable, Mapping


def open_text(path, mode="rt"):
    path = Path(path)
    return gzip.open(path, mode, encoding="utf-8", newline="") if path.suffix == ".gz" else open(path, mode, encoding="utf-8", newline="")


def read_tsv(path):
    with open_text(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows: Iterable[Mapping], fields=None):
    rows = list(rows)
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if fields is None:
        fields = []
        for row in rows:
            for key in row:
                if key not in fields:
                    fields.append(key)
    with open_text(path, "wt") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def read_fasta(path):
    sequences = {}
    key = None
    with open_text(path) as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                key = line[1:].split()[0]
                if key in sequences:
                    raise ValueError(f"Duplicate FASTA identifier: {key}")
                sequences[key] = ""
            elif key is None:
                raise ValueError("FASTA sequence encountered before a header")
            else:
                sequences[key] += line.upper().replace("T", "U")
    return sequences


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()
