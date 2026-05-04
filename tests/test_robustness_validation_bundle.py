import importlib.util
import sys
from pathlib import Path

from scipy.stats import hypergeom


SCRIPT_PATH = Path(__file__).resolve().parents[1] / "scripts" / "analysis" / "robustness_validation_bundle.py"
SPEC = importlib.util.spec_from_file_location("robustness_validation_bundle", SCRIPT_PATH)
p0 = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = p0
SPEC.loader.exec_module(p0)


def test_bh_adjust_preserves_missing_and_controls_monotonicity():
    adjusted = p0.bh_adjust([0.01, 0.04, None, 0.03])

    assert adjusted[2] is None
    assert adjusted[0] == 0.03
    assert adjusted[1] == 0.04
    assert adjusted[3] == 0.04


def test_hypergeom_support_test_counts_overlap_and_expected_value():
    result = p0.hypergeom_support_test(
        target_genes=["A", "B", "C"],
        validation_genes=["A", "C", "D", "E"],
        background_genes=list("ABCDEFGHIJ"),
    )

    assert result["background_n"] == 10
    assert result["validation_n_in_background"] == 4
    assert result["target_n_in_background"] == 3
    assert result["overlap_n"] == 2
    assert result["expected_overlap"] == 1.2
    assert result["overlap_genes"] == "A; C"
    assert result["p_value"] == hypergeom.sf(1, 10, 4, 3)


def test_build_gene_pubmed_query_includes_symbol_description_and_context():
    query = p0.build_gene_pubmed_query(
        gene="STAT3",
        description="Signal Transducer And Activator Of Transcription 3",
        context_term="oxidative stress",
    )

    assert '"STAT3"[Title/Abstract]' in query
    assert '"STAT3 gene"[Title/Abstract]' in query
    assert '"Signal Transducer And Activator Of Transcription 3"[Title/Abstract]' in query
    assert '"oxidative stress"[Title/Abstract]' in query


def test_build_gene_pubmed_query_avoids_bare_short_symbols():
    query = p0.build_gene_pubmed_query(
        gene="CA1",
        description="Carbonic Anhydrase 1",
        context_term="oxidative stress",
    )

    assert '"CA1"[Title/Abstract]' not in query
    assert '"CA1 gene"[Title/Abstract]' in query
    assert '"CA1 protein"[Title/Abstract]' in query
    assert '"Carbonic Anhydrase 1"[Title/Abstract]' in query
