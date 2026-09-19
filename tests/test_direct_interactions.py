from pathlib import Path

from importlib.util import module_from_spec, spec_from_file_location


SCRIPT = Path(__file__).parents[1] / "benchmarks" / "direct_interactions" / "prepare_direct_interactions.py"
SPEC = spec_from_file_location("prepare_direct_interactions", SCRIPT)
MODULE = module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_transcript_region_boundaries():
    model = ("coding", 100, 400)
    assert MODULE.transcript_region(model, 10, 20) == "5'UTR"
    assert MODULE.transcript_region(model, 200, 220) == "CDS"
    assert MODULE.transcript_region(model, 450, 470) == "3'UTR"


def test_transcript_region_unmapped_and_noncoding():
    assert MODULE.transcript_region(None, 10, 20) == "mRNA_region_unmapped"
    assert MODULE.transcript_region(("noncoding_exonic", None, None), 10, 20) == "noncoding_exonic"
