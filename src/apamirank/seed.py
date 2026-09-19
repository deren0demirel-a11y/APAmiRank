from __future__ import annotations


RNA_COMPLEMENT = str.maketrans("ACGUTN", "UGCAAN")


def normalize_mirna(sequence: str) -> str:
    seq = "".join(sequence.split()).upper().replace("T", "U")
    invalid = set(seq) - set("ACGUN")
    if invalid:
        raise ValueError(f"Invalid miRNA characters: {sorted(invalid)}")
    if len(seq) < 8:
        raise ValueError("miRNA sequence must contain at least 8 nucleotides")
    if "N" in seq[:8]:
        raise ValueError("miRNA nucleotides 1-8 cannot contain N")
    return seq


def reverse_complement_rna(sequence: str) -> str:
    return sequence.translate(RNA_COMPLEMENT)[::-1]


def canonical_target_motifs(mirna_sequence: str):
    """Return canonical target motifs in target-RNA 5'->3' orientation.

    miRNA positions use the conventional 1-based definition. The target A1
    adenosine is appended at the 3' end of the reverse-complemented seed when
    written in target 5'->3' orientation.
    """
    seq = normalize_mirna(mirna_sequence)
    seed_2_8 = reverse_complement_rna(seq[1:8])
    seed_2_7 = reverse_complement_rna(seq[1:7])
    return [
        ("8mer", seed_2_8 + "A", 4),
        ("7mer-m8", seed_2_8, 3),
        ("7mer-A1", seed_2_7 + "A", 2),
        ("6mer", seed_2_7, 1),
    ]
