#!/usr/bin/env python3
"""Unit tests for PubMedConnector module."""

import json
import sys
import unittest
from pathlib import Path
from typing import Set

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src.connectors.pubmed_connector import (
    PubMedConnector,
    PubMedConfig,
    PubMedError,
    build_cooccurrence_query,
    build_search_term,
)


class TestBuildSearchTerm(unittest.TestCase):
    """Test search term construction."""

    def test_preserve_full_pathway_name(self):
        """Full pathway name should be preserved as exact phrase."""
        q = build_search_term("MAPK signaling pathway")
        self.assertEqual(q, '"MAPK signaling pathway"[Title/Abstract]')

    def test_no_blind_strip_of_pathway_word(self):
        """Old bug: stripping 'pathway' broke exact-phrase matching."""
        q = build_search_term("C-type lectin receptor signaling pathway")
        self.assertIn("C-type lectin receptor signaling pathway", q)
        # The word 'pathway' must remain in the query (not stripped)
        self.assertIn("pathway", q)

    def test_age_rage_expansion(self):
        """AGE-RAGE acronym should expand to full name."""
        q = build_search_term("AGE-RAGE signaling pathway in diabetic complications")
        self.assertIn("advanced glycation end product", q)
        self.assertIn("AGE-RAGE", q)

    def test_quote_escaping(self):
        """Quotes inside pathway names should be stripped."""
        q = build_search_term('MAPK "special" pathway')
        self.assertNotIn('"special"', q)  # inner quotes removed


class TestBuildCooccurrenceQuery(unittest.TestCase):
    """Test co-occurrence query construction."""

    def test_basic_cooccurrence(self):
        q = build_cooccurrence_query("MAPK signaling pathway", "oxidative stress")
        self.assertIn("MAPK signaling pathway", q)
        self.assertIn("oxidative stress", q)
        self.assertIn("AND", q)

    def test_context_term_escaping(self):
        q = build_cooccurrence_query("X", 'oxidative "stress"')
        self.assertNotIn('"stress"', q)
        self.assertIn("oxidative stress", q)


class TestPubMedConnectorCache(unittest.TestCase):
    """Test caching behaviour."""

    def setUp(self):
        self.conn = PubMedConnector(cache_dir=Path(__file__).parent / "_test_cache")
        self.conn.clear_cache()

    def tearDown(self):
        self.conn.clear_cache()

    def test_cache_starts_empty(self):
        self.assertEqual(len(self.conn._cache), 0)

    def test_manual_cache_seed(self):
        self.conn._cache["test::oxidative stress"] = {
            "count": 42,
            "timestamp": "2026-04-26",
            "query": "test",
        }
        self.conn._save_cache()
        # Fresh connector should load the cache
        conn2 = PubMedConnector(cache_dir=Path(__file__).parent / "_test_cache")
        self.assertEqual(len(conn2._cache), 1)
        self.assertEqual(conn2._cache["test::oxidative stress"]["count"], 42)


class TestComputeMetrics(unittest.TestCase):
    """Test offline metric computation against known validation data."""

    @classmethod
    def setUpClass(cls):
        validation_path = (
            Path(__file__).resolve().parent.parent
            / "src" / "outputs" / "sensitivity" / "pubmed_validation.json"
        )
        if not validation_path.exists():
            raise unittest.SkipTest("pubmed_validation.json not found")

        with open(validation_path, "r", encoding="utf-8") as f:
            cls.data = json.load(f)

    def test_metrics_match_known_results(self):
        """Computed metrics must match the previously reported values."""
        details = self.data["pathway_details"]
        pub_pos: Set[str] = {d["pathway"] for d in details if d["pubmed_count"] >= 10}
        fw_pos: Set[str] = {d["pathway"] for d in details if d["framework_positive"]}
        all_p: Set[str] = {d["pathway"] for d in details}

        m = PubMedConnector.compute_metrics(pub_pos, fw_pos, all_p)

        self.assertAlmostEqual(m["precision"], 0.8438, places=4)
        self.assertAlmostEqual(m["recall"], 0.9643, places=4)
        self.assertAlmostEqual(m["f1"], 0.9000, places=4)
        self.assertAlmostEqual(m["jaccard"], 0.8182, places=4)

    def test_batch_cooccurrence_with_seeded_cache(self):
        """Batch query should return correct counts when cache is seeded."""
        conn = PubMedConnector(cache_dir=Path(__file__).parent / "_test_cache")
        conn.clear_cache()

        # Seed from validation JSON
        for d in self.data["pathway_details"]:
            p = d["pathway"]
            c = d["pubmed_count"]
            conn._cache[f"{p}::oxidative stress"] = {
                "count": c,
                "timestamp": "2026-04-26",
                "query": f'"{p}"[Title/Abstract] AND "oxidative stress"[Title/Abstract]',
            }
        conn._save_cache()

        counts = conn.batch_cooccurrence(
            ["MAPK signaling pathway",
             "AGE-RAGE signaling pathway in diabetic complications",
             "Nitrogen metabolism"],
            show_progress=False,
        )
        self.assertEqual(counts["MAPK signaling pathway"], 7091)
        self.assertEqual(counts["AGE-RAGE signaling pathway in diabetic complications"], 0)
        self.assertEqual(counts["Nitrogen metabolism"], 320)

        conn.clear_cache()


if __name__ == "__main__":
    unittest.main(verbosity=2)
