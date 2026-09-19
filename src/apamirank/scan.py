from __future__ import annotations

import re
from collections import defaultdict

from .io import read_fasta, read_tsv, write_tsv
from .seed import canonical_target_motifs, normalize_mirna


REQUIRED_METADATA = {"isoform_id", "gene", "locus_id", "pas_rank", "chrom", "start0", "end", "strand"}


def _as_int(row, field):
    try:
        return int(row[field])
    except Exception as exc:
        raise ValueError(f"Invalid integer in {field}: {row.get(field)!r}") from exc


def scan_sites(metadata_path, fasta_path, mirna_id, mirna_sequence, output_path, eligible_genes=None, flank=50):
    meta = read_tsv(metadata_path)
    if not meta:
        raise ValueError("Isoform metadata is empty")
    missing = REQUIRED_METADATA - set(meta[0])
    if missing:
        raise ValueError("Missing metadata columns: " + ", ".join(sorted(missing)))
    sequences = read_fasta(fasta_path)
    allowed = None
    if eligible_genes:
        with open(eligible_genes, encoding="utf-8") as handle:
            allowed = {line.strip() for line in handle if line.strip() and not line.startswith("#")}

    loci = defaultdict(list)
    seen_isoforms = set()
    for row in meta:
        iso = row["isoform_id"]
        if iso in seen_isoforms:
            raise ValueError(f"Duplicate isoform_id in metadata: {iso}")
        seen_isoforms.add(iso)
        if iso not in sequences:
            raise ValueError(f"Metadata isoform missing from FASTA: {iso}")
        if allowed is not None and row["gene"] not in allowed:
            continue
        start, end = _as_int(row, "start0"), _as_int(row, "end")
        if end <= start:
            raise ValueError(f"Invalid interval for {iso}: {start}-{end}")
        if len(sequences[iso]) != end - start:
            raise ValueError(f"Sequence length does not match genomic interval for {iso}")
        if row["strand"] not in {"+", "-"}:
            raise ValueError(f"Invalid strand for {iso}: {row['strand']}")
        row = dict(row)
        row["pas_rank"] = _as_int(row, "pas_rank")
        row["start0"] = start
        row["end"] = end
        row["sequence"] = sequences[iso]
        loci[(row["gene"], row["locus_id"])].append(row)

    motifs = canonical_target_motifs(mirna_sequence)
    mirna_sequence = normalize_mirna(mirna_sequence)
    output = []
    for (gene, locus_id), isoforms in sorted(loci.items()):
        ranks = sorted(r["pas_rank"] for r in isoforms)
        if ranks != list(range(1, len(ranks) + 1)):
            raise ValueError(f"PAS ranks must be unique and contiguous within {gene}/{locus_id}; observed {ranks}")
        hits = {}
        for row in sorted(isoforms, key=lambda r: r["pas_rank"]):
            seq = row["sequence"]
            occupied = []
            for seed_class, motif, strength in motifs:
                for match in re.finditer(f"(?={motif})", seq):
                    start = match.start()
                    stop = start + len(motif)
                    if any(start >= left and stop <= right for left, right in occupied):
                        continue
                    occupied.append((start, stop))
                    if row["strand"] == "+":
                        genomic_start, genomic_end = row["start0"] + start, row["start0"] + stop
                    else:
                        genomic_start, genomic_end = row["end"] - stop, row["end"] - start
                    key = (row["chrom"], genomic_start, genomic_end, row["strand"], seed_class)
                    entry = hits.setdefault(key, {
                        "ranks": set(), "isoforms": set(), "contexts": [],
                        "strength": strength, "motif": motif,
                    })
                    entry["ranks"].add(row["pas_rank"])
                    entry["isoforms"].add(row["isoform_id"])
                    entry["contexts"].append((row["pas_rank"], seq, start))

        n_isoforms = len(isoforms)
        for (chrom, genomic_start, genomic_end, strand, seed_class), hit in hits.items():
            present = sorted(hit["ranks"])
            first, last = min(present), max(present)
            # Use the longest isoform carrying the site so downstream AU context
            # is not artificially truncated at a proximal PAS.
            _, seq, position = max(hit["contexts"], key=lambda item: item[0])
            motif = hit["motif"]
            window_start = max(0, position - flank)
            window_end = min(len(seq), position + len(motif) + flank)
            flanking = seq[max(0, position - 30):position] + seq[position + len(motif):min(len(seq), position + len(motif) + 30)]
            au = sum(base in "AU" for base in flanking) / len(flanking) if flanking else 0.0
            output.append({
                "miRNA_id": mirna_id,
                "miRNA_sequence_5to3": mirna_sequence,
                "gene": gene,
                "locus_id": locus_id,
                "seed_class": seed_class,
                "site_sequence_target_RNA_5to3": motif,
                "chrom": chrom,
                "genomic_start0": genomic_start,
                "genomic_end": genomic_end,
                "strand": strand,
                "genomic_coordinate_1based": f"{chrom}:{genomic_start + 1}-{genomic_end}",
                "first_available_PAS_rank": first,
                "last_available_PAS_rank": last,
                "present_PAS_ranks": ";".join(map(str, present)),
                "n_isoforms_present": len(present),
                "n_isoforms_locus": n_isoforms,
                "region_class": "common" if first == 1 else ("distal" if first == n_isoforms else "intermediate"),
                "extension_specific_yes_no": "No" if first == 1 else "Yes",
                "retained_in_all_longer_isoforms": "Yes" if present == list(range(first, n_isoforms + 1)) else "No",
                "seed_strength_ordinal": hit["strength"],
                "AU_flank30": f"{au:.6f}",
                "window_start1_UTR": window_start + 1,
                "window_end1_UTR": window_end,
                "site_start1_window": position - window_start + 1,
                "site_end1_window": position - window_start + len(motif),
                "window_sequence_RNA_5to3": seq[window_start:window_end],
                "source_isoforms": ";".join(sorted(hit["isoforms"])),
            })

    output.sort(key=lambda r: (r["gene"], r["locus_id"], r["chrom"], int(r["genomic_start0"]), -int(r["seed_strength_ordinal"])))
    safe_id = re.sub(r"[^A-Za-z0-9]+", "_", mirna_id).strip("_").upper() or "MIRNA"
    for index, row in enumerate(output, 1):
        row["site_id"] = f"APAMIRANK_{safe_id}_{index:08d}"
    fields = ["site_id"] + [field for field in output[0] if field != "site_id"] if output else ["site_id"]
    write_tsv(output_path, output, fields)
    return output
