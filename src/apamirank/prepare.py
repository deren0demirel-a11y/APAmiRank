from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

from .io import write_tsv


DNA_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def prepare_qapa(qapa_bed, genome_fasta, metadata_out, fasta_out):
    try:
        from pyfaidx import Fasta
    except ImportError as exc:
        raise RuntimeError("pyfaidx is required for prepare-qapa") from exc
    genome = Fasta(genome_fasta, as_raw=True, sequence_always_upper=True)
    raw = []
    with open(qapa_bed, encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 7:
                raise ValueError(f"QAPA BED line {line_number} has fewer than 7 columns")
            chrom, start, end, name, _, strand, gene = parts[:7]
            start, end = int(start), int(end)
            locus_key = (gene, chrom, strand, start if strand == "+" else end)
            raw.append({"chrom": chrom, "start0": start, "end": end, "name": name, "strand": strand, "gene": gene, "locus_key": locus_key})
    groups = defaultdict(list)
    for row in raw:
        groups[row["locus_key"]].append(row)
    metadata, sequences = [], []
    used = set()
    for locus_number, (key, rows) in enumerate(sorted(groups.items()), 1):
        gene, chrom, strand, _ = key
        locus_id = f"{gene}_L{locus_number}"
        rows.sort(key=lambda row: (row["end"] - row["start0"], row["start0"], row["end"], row["name"]))
        for rank, row in enumerate(rows, 1):
            isoform_id = row["name"] or f"{locus_id}_I{rank}"
            if isoform_id in used:
                isoform_id = f"{isoform_id}__{locus_id}_I{rank}"
            used.add(isoform_id)
            try:
                sequence = str(genome[row["chrom"]][row["start0"]:row["end"]])
            except Exception as exc:
                raise ValueError(f"Cannot extract {row['chrom']}:{row['start0']}-{row['end']}") from exc
            if strand == "-":
                sequence = sequence.translate(DNA_COMPLEMENT)[::-1]
            metadata.append({
                "isoform_id": isoform_id, "gene": gene, "locus_id": locus_id,
                "pas_rank": rank, "chrom": chrom, "start0": row["start0"],
                "end": row["end"], "strand": strand,
            })
            sequences.append((isoform_id, sequence))
    write_tsv(metadata_out, metadata)
    fasta_out = Path(fasta_out)
    fasta_out.parent.mkdir(parents=True, exist_ok=True)
    with open(fasta_out, "w", encoding="utf-8") as handle:
        for isoform_id, sequence in sequences:
            handle.write(f">{isoform_id}\n{sequence}\n")
    return metadata
