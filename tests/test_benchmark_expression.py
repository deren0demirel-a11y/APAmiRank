import csv
import gzip
import importlib.util
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / "benchmarks" / "gse52531" / "prepare_expression.py"
SPEC = importlib.util.spec_from_file_location("prepare_expression", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


def test_expression_converter_preserves_replicates(tmp_path):
    source = tmp_path / "input.txt.gz"
    output = tmp_path / "output.tsv"
    with gzip.open(source, "wt", encoding="utf-8") as handle:
        handle.write(
            "#geneID\tgeneLen\tPuc19_124_1,Puc19_124_2\t"
            "Puc19_155_1,Puc19_155_2\tmiR-124_1,miR-124_2\tmiR-155_1,miR-155_2\n"
            "NM_TEST\t1000\t1.0,2.0\t3.0,4.0\t0.5,0.6\t7.0,8.0\n"
        )

    count = MODULE.convert(source, output, "Huh7")

    assert count == 1
    with output.open(encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 8
    assert {row["condition"] for row in rows} == {
        "mock_miR-124",
        "mock_miR-155",
        "miR-124",
        "miR-155",
    }
    assert rows[0]["sample_id"] == "Puc19_124_1"
    assert rows[-1]["rpkm"] == "8.0"
