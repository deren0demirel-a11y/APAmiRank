import csv
import importlib.util
import json
import sys
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / "benchmarks" / "gse52531" / "prepare_3p_usage.py"
SPEC = importlib.util.spec_from_file_location("prepare_3p_usage", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def test_conservative_strand_aware_pas_mapping(tmp_path):
    refgene = tmp_path / "refGene.txt"
    bed = tmp_path / "clusters.bed"
    qapa = tmp_path / "qapa.bed"
    usage = tmp_path / "usage.tsv"
    mapping = tmp_path / "map.tsv"
    qc_path = tmp_path / "qc.json"
    refgene.write_text(
        "0\tNM_PLUS\tchr1\t+\t100\t300\t120\t200\t2\t100,180,\t150,300,\t0\tPLUS\n"
        "0\tNM_MINUS\tchr2\t-\t500\t800\t600\t760\t2\t500,700,\t650,800,\t0\tMINUS\n",
        encoding="utf-8",
    )
    bed.write_text(
        "chr1\t205\t220\t3\t0\t+\n"
        "chr1\t255\t270\t1\t0\t+\n"
        "chr1\t275\t290\t7\t0\t+\n"
        "chr2\t580\t595\t4\t0\t-\n"
        "chr2\t520\t535\t6\t0\t-\n",
        encoding="utf-8",
    )

    loci = MODULE.load_refgene(refgene)
    counts, lookup, qc = MODULE.map_clusters(bed, loci, min_count=2)
    report = MODULE.write_outputs(
        counts, lookup, "test", qapa, usage, mapping, qc_path, qc
    )

    assert report["retained_apa_loci"] == 2
    assert report["retained_isoforms"] == 4
    assert report["min_locus_count"] == 10.0
    assert report["below_min_count"] == 1
    with usage.open(encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    plus = [row for row in rows if row["isoform_id"].startswith("PLUS__")]
    assert [row["isoform_id"].rsplit("PAS", 1)[1] for row in plus] == ["220", "290"]
    assert abs(sum(float(row["usage_fraction"]) for row in plus) - 1.0) < 1e-10
    assert json.loads(qc_path.read_text())["pas_definition"].startswith("BED end")
