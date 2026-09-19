from pathlib import Path

from apamirank.scan import scan_sites


ROOT = Path(__file__).resolve().parents[1]


def test_scanner_detects_common_and_extension_specific_sites(tmp_path):
    output = tmp_path / "sites.tsv"
    rows = scan_sites(
        ROOT / "examples/minimal/isoforms.tsv",
        ROOT / "examples/minimal/isoforms.fa",
        "hsa-miR-example",
        "AUGCAUGCAUGCAUGCAUGCAU",
        output,
    )
    assert output.exists()
    assert rows
    assert {row["gene"] for row in rows} == {"GENEA", "GENEB"}
    gene_a = [row for row in rows if row["gene"] == "GENEA"]
    assert any(row["region_class"] == "common" for row in gene_a)
    assert any(row["extension_specific_yes_no"] == "Yes" for row in gene_a)
    assert all(row["miRNA_id"] == "hsa-miR-example" for row in rows)
