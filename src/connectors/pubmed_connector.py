"""
PubMed / NCBI E-utilities Connector Module

RESTful API for PubMed literature search:
https://eutils.ncbi.nlm.nih.gov/entrez/eutils/

Supports: article count queries, co-occurrence validation, local JSON caching.
"""

import json
import logging
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Set

import requests

logger = logging.getLogger(__name__)

PUBMED_ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"


# ---------------------------------------------------------------------------
# Known search-term expansions for acronyms that fail when stripped
# ---------------------------------------------------------------------------
_SEARCH_EXPANSIONS: Dict[str, str] = {
    "age-rage": '"advanced glycation end product" OR "AGE-RAGE"',
}


def build_search_term(pathway_name: str) -> str:
    """
    Build a PubMed Title/Abstract search term from a pathway name.

    Strategy:
    1.  Keep the full pathway name as the primary phrase (do NOT blindly
        strip "pathway" / "signaling" — that destroys exact-phrase matching).
    2.  For known acronyms (e.g. AGE-RAGE), expand to the full name.
    3.  Escape quotes to avoid query syntax errors.

    Examples
    --------
    >>> build_search_term("MAPK signaling pathway")
    '"MAPK signaling pathway"[Title/Abstract]'

    >>> build_search_term("AGE-RAGE signaling pathway in diabetic complications")
    '("advanced glycation end product" OR "AGE-RAGE")[Title/Abstract]'
    """
    lower = pathway_name.lower()

    # Check for known acronym expansions first
    for key, expansion in _SEARCH_EXPANSIONS.items():
        if key in lower:
            return f"({expansion})[Title/Abstract]"

    # Default: use the full pathway name as an exact phrase
    safe = pathway_name.replace('"', "")
    return f'"{safe}"[Title/Abstract]'


def build_cooccurrence_query(pathway_name: str, context_term: str = "oxidative stress") -> str:
    """
    Build a PubMed query for co-occurrence of *pathway_name* and *context_term*
    in Title/Abstract.
    """
    pathway_term = build_search_term(pathway_name)
    safe_context = context_term.replace('"', "")
    return f"({pathway_term}) AND (\"{safe_context}\"[Title/Abstract])"


@dataclass(frozen=True)
class PubMedConfig:
    """PubMed E-utilities configuration."""

    base_url: str = PUBMED_ESEARCH_URL
    timeout: int = 30
    max_retries: int = 3
    delay_seconds: float = 0.35  # NCBI recommends ~3 req/sec
    user_agent: str = "Mozilla/5.0 (compatible; ResearchBot/1.0)"


class PubMedError(Exception):
    """PubMed API error."""
    pass


class PubMedConnector:
    """
    NCBI E-utilities client with local JSON caching.

    Cache file: ``<cache_dir>/pubmed_cache.json``
    Keys: ``"<pathway_name>::<context_term>"``
    Values: ``{"count": int, "timestamp": str}``
    """

    def __init__(
        self,
        config: Optional[PubMedConfig] = None,
        cache_dir: Optional[Path] = None,
    ):
        self.config = config or PubMedConfig()
        self._session = requests.Session()
        self._session.headers.update({"User-Agent": self.config.user_agent})

        # Cache setup
        if cache_dir is None:
            cache_dir = Path(__file__).resolve().parent.parent / "outputs" / "sensitivity"
        self._cache_dir = Path(cache_dir)
        self._cache_dir.mkdir(parents=True, exist_ok=True)
        self._cache_path = self._cache_dir / "pubmed_cache.json"
        self._cache: Dict[str, Dict] = self._load_cache()

    # ------------------------------------------------------------------
    # Cache helpers
    # ------------------------------------------------------------------

    def _load_cache(self) -> Dict[str, Dict]:
        if self._cache_path.exists():
            try:
                with open(self._cache_path, "r", encoding="utf-8") as f:
                    return json.load(f)
            except (json.JSONDecodeError, OSError) as e:
                logger.warning(f"PubMed cache corrupt, resetting: {e}")
        return {}

    def _save_cache(self) -> None:
        with open(self._cache_path, "w", encoding="utf-8") as f:
            json.dump(self._cache, f, indent=2, ensure_ascii=False)

    def _cache_key(self, pathway_name: str, context_term: str) -> str:
        return f"{pathway_name}::{context_term}"

    def clear_cache(self) -> None:
        """Remove all cached entries."""
        self._cache.clear()
        if self._cache_path.exists():
            self._cache_path.unlink()
        logger.info("PubMed cache cleared.")

    # ------------------------------------------------------------------
    # Core API
    # ------------------------------------------------------------------

    def search_count(self, query: str, use_cache: bool = True) -> int:
        """
        Return the number of PubMed articles matching *query*.

        Parameters
        ----------
        query : str
            Raw PubMed query string.
        use_cache : bool
            If True, return cached result when available.

        Returns
        -------
        int
            Article count (>= 0).  Returns -1 on API failure.
        """
        cache_key = query
        if use_cache and cache_key in self._cache:
            logger.debug(f"Cache hit: {query[:60]}...")
            return self._cache[cache_key]["count"]

        params = {
            "db": "pubmed",
            "term": query,
            "retmax": 0,
            "retmode": "json",
        }

        last_error: Optional[Exception] = None
        for attempt in range(1, self.config.max_retries + 1):
            try:
                resp = self._session.get(
                    self.config.base_url,
                    params=params,
                    timeout=self.config.timeout,
                )
                resp.raise_for_status()
                data = resp.json()
                count = int(data.get("esearchresult", {}).get("count", 0))

                # Cache and rate-limit
                if use_cache:
                    self._cache[cache_key] = {
                        "count": count,
                        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S"),
                        "query": query,
                    }
                    self._save_cache()

                time.sleep(self.config.delay_seconds)
                return count

            except requests.RequestException as e:
                last_error = e
                logger.warning(
                    f"PubMed request attempt {attempt}/{self.config.max_retries} failed: {e}"
                )
                time.sleep(self.config.delay_seconds * attempt)

        logger.error(f"PubMed query failed after all retries: {last_error}")
        return -1

    def cooccurrence_count(
        self,
        pathway_name: str,
        context_term: str = "oxidative stress",
        use_cache: bool = True,
    ) -> int:
        """
        Count PubMed articles where *pathway_name* co-occurs with
        *context_term* in Title/Abstract.

        Parameters
        ----------
        pathway_name : str
            KEGG / GO pathway name (e.g. "MAPK signaling pathway").
        context_term : str
            Context keyword for co-occurrence (default: "oxidative stress").
        use_cache : bool

        Returns
        -------
        int
            Co-occurrence article count, or -1 on API failure.
        """
        cache_key = self._cache_key(pathway_name, context_term)
        if use_cache and cache_key in self._cache:
            logger.debug(f"Cache hit: {pathway_name}")
            return self._cache[cache_key]["count"]

        query = build_cooccurrence_query(pathway_name, context_term)
        count = self.search_count(query, use_cache=False)  # inner caching handled below

        if count >= 0 and use_cache:
            self._cache[cache_key] = {
                "count": count,
                "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S"),
                "pathway_name": pathway_name,
                "context_term": context_term,
                "query": query,
            }
            self._save_cache()

        return count

    def batch_cooccurrence(
        self,
        pathway_names: List[str],
        context_term: str = "oxidative stress",
        use_cache: bool = True,
        show_progress: bool = True,
    ) -> Dict[str, int]:
        """
        Query co-occurrence counts for multiple pathways.

        Returns
        -------
        Dict[str, int]
            Mapping ``pathway_name -> count``.  Failed queries return -1.
        """
        results: Dict[str, int] = {}
        total = len(pathway_names)

        for i, name in enumerate(pathway_names, 1):
            count = self.cooccurrence_count(name, context_term, use_cache=use_cache)
            results[name] = count

            if show_progress:
                status = "[+]" if count >= 10 else "[ ]"
                print(f"  [{i:3d}/{total}] {status} {name[:45]:<45} n={count:4d}")

        return results

    # ------------------------------------------------------------------
    # Validation helpers
    # ------------------------------------------------------------------

    @staticmethod
    def compute_metrics(
        pubmed_positive: Set[str],
        framework_positive: Set[str],
        all_pathways: Set[str],
    ) -> Dict[str, float]:
        """
        Compute Precision, Recall, F1, Jaccard, Accuracy.

        Parameters
        ----------
        pubmed_positive : Set[str]
            Pathways deemed positive by PubMed gold standard.
        framework_positive : Set[str]
            Pathways deemed positive by the classification framework.
        all_pathways : Set[str]
            Universe of all tested pathways.

        Returns
        -------
        Dict[str, float]
            {'tp', 'fp', 'fn', 'tn', 'precision', 'recall', 'f1', 'jaccard', 'accuracy'}
        """
        tp = len(pubmed_positive & framework_positive)
        fp = len(framework_positive - pubmed_positive)
        fn = len(pubmed_positive - framework_positive)
        tn = len(all_pathways - pubmed_positive - framework_positive)

        precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
        recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
        f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0.0
        jaccard = tp / (tp + fp + fn) if (tp + fp + fn) > 0 else 0.0
        accuracy = (tp + tn) / len(all_pathways) if all_pathways else 0.0

        return {
            "tp": tp,
            "fp": fp,
            "fn": fn,
            "tn": tn,
            "precision": precision,
            "recall": recall,
            "f1": f1,
            "jaccard": jaccard,
            "accuracy": accuracy,
        }
