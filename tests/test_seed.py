from apamirank.seed import canonical_target_motifs, normalize_mirna


def test_normalization_and_generic_seed_construction():
    assert normalize_mirna("augcaugcaugc") == "AUGCAUGCAUGC"
    motifs = {name: motif for name, motif, _ in canonical_target_motifs("AUGCAUGCAUGC")}
    assert motifs == {
        "8mer": "GCAUGCAA",
        "7mer-m8": "GCAUGCA",
        "7mer-A1": "CAUGCAA",
        "6mer": "CAUGCA",
    }


def test_different_mirnas_produce_different_motifs():
    first = canonical_target_motifs("AUGCAUGCAUGC")
    second = canonical_target_motifs("UGAGGUAGUAGGUUGUAUAGUU")
    assert first != second
