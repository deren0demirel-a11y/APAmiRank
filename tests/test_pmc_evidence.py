import importlib.util
import xml.etree.ElementTree as ET
from pathlib import Path


MODULE_PATH = (
    Path(__file__).parents[1]
    / "benchmarks"
    / "locked_luciferase_holdout"
    / "extract_pmc_evidence.py"
)
SPEC = importlib.util.spec_from_file_location("extract_pmc_evidence", MODULE_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


def test_evidence_units_requires_reporter_and_design_language():
    article = ET.fromstring(
        """<article><body><sec><title>Results</title>
        <p>The wild-type 3′UTR luciferase reporter was repressed, whereas the
        mutant construct was not.</p>
        <p>Western blot showed lower protein abundance.</p>
        </sec></body></article>"""
    )
    units = MODULE.evidence_units(article)
    assert units == [
        (
            "Results",
            "The wild-type 3′UTR luciferase reporter was repressed, whereas the mutant construct was not.",
        )
    ]


def test_article_id_normalizes_pmc_prefix_later_in_main():
    article = ET.fromstring(
        '<article><article-meta><article-id pub-id-type="pmc">123</article-id>'
        '<article-id pub-id-type="pmid">456</article-id></article-meta></article>'
    )
    assert MODULE.article_id(article, "pmc") == "123"
    assert MODULE.article_id(article, "pmid") == "456"
