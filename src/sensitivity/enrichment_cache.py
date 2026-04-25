"""
富集分析缓存工具

避免对同一基因集反复调用富集API，通过MD5(genes|pvalue)作为缓存key，
将结果序列化为JSON存至 `src/outputs/sensitivity/enrichment_cache/`。
"""

from pathlib import Path
import hashlib
import json
from typing import Dict, Optional

import pandas as pd

from src.config import BACKGROUND_GENES
from src.enrichment.enrichment_analysis import run_enrichment_analysis

# ---------------------------------------------------------------------------
# 缓存目录
# ---------------------------------------------------------------------------
CACHE_DIR = Path(__file__).resolve().parent.parent.parent / "src" / "outputs" / "sensitivity" / "enrichment_cache"
CACHE_DIR.mkdir(parents=True, exist_ok=True)

# 默认参数
PVALUE_CUTOFF = 0.01
TOP_N = 20


def _cache_key(gene_list, pvalue_cutoff):
    genes_str = ";".join(sorted(gene_list))
    return hashlib.md5(f"{genes_str}|{pvalue_cutoff}".encode()).hexdigest()[:12]


def _load(key: str) -> Optional[Dict]:
    f = CACHE_DIR / f"{key}.json"
    if not f.exists():
        return None
    with open(f, "r") as fh:
        data = json.load(fh)
    return {k: pd.DataFrame(v) if isinstance(v, list) else v for k, v in data.items()}


def _save(key: str, significant: Dict):
    f = CACHE_DIR / f"{key}.json"
    serializable = {
        k: (v.to_dict(orient="records") if isinstance(v, pd.DataFrame) else v)
        for k, v in significant.items()
    }
    with open(f, "w") as fh:
        json.dump(serializable, fh)


def cached_enrichment(gene_list, pvalue_cutoff=PVALUE_CUTOFF, background=BACKGROUND_GENES, top_n=TOP_N):
    """带缓存的富集分析。

    相同gene_list + pvalue_cutoff组合只会调用一次API，
    后续直接从本地JSON读取。
    """
    key = _cache_key(gene_list, pvalue_cutoff)
    cached = _load(key)
    if cached is not None:
        return {"status": "成功", "significant": cached}

    enrichment = run_enrichment_analysis(
        gene_list=gene_list,
        pvalue_cutoff=pvalue_cutoff,
        top_n=top_n,
        background=background,
    )
    if enrichment["status"] == "成功" and "significant" in enrichment:
        _save(key, enrichment["significant"])
    return enrichment
