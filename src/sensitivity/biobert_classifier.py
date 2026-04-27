"""
BioBERT / Sentence-BERT 语义通路分类器

利用 sentence-transformers 的预训练语义模型，
将通路名称与疾病概念编码为稠密向量，通过余弦相似度判定语义相关性。

模型默认使用 all-MiniLM-L6-v2（~80MB，专门训练用于语义相似度）。
如需更强的生物医学语义理解，可替换为生物医学专用模型。

特性
----
- 零样本：无需标注数据、无需微调
- 懒加载：首次调用时才加载模型
- 缓存：编码结果持久化到 pickle，二次运行秒开
- 批量编码：自动分批处理，适配 CPU/GPU
"""

import hashlib
import os
import pickle
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import numpy as np

# 阻止 transformers 尝试导入有问题的 TensorFlow
os.environ.setdefault("USE_TF", "0")
os.environ.setdefault("TRANSFORMERS_NO_TF", "1")

# ---------------------------------------------------------------------------
# 常量
# ---------------------------------------------------------------------------
# 默认模型：轻量、专门为语义相似度优化的通用模型
DEFAULT_MODEL = "sentence-transformers/all-MiniLM-L6-v2"

# 生物医学专用备选（更大、更专业，需约 400MB）
BIOMEDICAL_MODEL = "pritamdeka/S-PubMedBert-MS-MARCO"

REFERENCE_TEXT = "oxidative stress"
DEFAULT_THRESHOLD = 0.55  # sentence-transformer 模型的分数分布更合理
BATCH_SIZE = 64

# 缓存路径
DEFAULT_CACHE_PATH = (
    Path(__file__).resolve().parent.parent
    / "outputs"
    / "sensitivity"
    / "biobert_cache.pkl"
)


class BioBERTClassifier:
    """Sentence-BERT 语义相似度分类器。"""

    def __init__(
        self,
        model_name: str = DEFAULT_MODEL,
        reference_text: str = REFERENCE_TEXT,
        cache_path: Optional[Union[str, Path]] = None,
        device: Optional[str] = None,
    ):
        """
        参数
        ----------
        model_name : str
            sentence-transformers 模型名称。
        reference_text : str
            疾病参考概念（默认 "oxidative stress"）。
        cache_path : str | Path | None
            编码缓存文件路径。None 则使用默认路径。
        device : str | None
            "cuda" / "cpu" / None（自动检测）。
        """
        self.model_name = model_name
        self.reference_text = reference_text
        self.cache_path = Path(cache_path) if cache_path else DEFAULT_CACHE_PATH
        self.device = device or ("cuda" if self._has_cuda() else "cpu")

        self._model = None
        self._ref_vector: Optional[np.ndarray] = None
        self._cache: Dict[str, np.ndarray] = {}

        # 尝试加载已有缓存
        self._load_cache()

    @staticmethod
    def _has_cuda() -> bool:
        try:
            import torch
            return torch.cuda.is_available()
        except Exception:
            return False

    # -----------------------------------------------------------------------
    # 模型加载
    # -----------------------------------------------------------------------
    def _load_model(self):
        """懒加载 SentenceTransformer 模型（仅调用一次）。"""
        if self._model is not None:
            return

        from sentence_transformers import SentenceTransformer

        # 国内网络 fallback
        try:
            self._model = SentenceTransformer(self.model_name, device=self.device)
        except Exception:
            os.environ["HF_ENDPOINT"] = "https://hf-mirror.com"
            self._model = SentenceTransformer(self.model_name, device=self.device)

        # 预编码参考概念
        self._ref_vector = self._encode_single(self.reference_text)

    def _encode_single(self, text: str) -> np.ndarray:
        """编码单条文本，返回 L2 归一化后的向量。"""
        self._load_model()
        vec = self._model.encode(text, convert_to_numpy=True, normalize_embeddings=True)
        return vec.astype(np.float32)

    # -----------------------------------------------------------------------
    # 批量编码
    # -----------------------------------------------------------------------
    def encode(self, texts: List[str]) -> np.ndarray:
        """
        批量编码文本列表。

        返回
        -------
        np.ndarray, shape (N, D)
            L2 归一化后的向量矩阵。D 取决于模型（MiniLM 为 384 维）。
        """
        self._load_model()

        vectors = []
        uncached_texts = []
        uncached_indices = []

        # 先从缓存取
        for i, text in enumerate(texts):
            key = self._cache_key(text)
            if key in self._cache:
                vectors.append((i, self._cache[key]))
            else:
                uncached_texts.append(text)
                uncached_indices.append(i)

        # 批量编码未缓存项
        if uncached_texts:
            batch_vecs = self._model.encode(
                uncached_texts,
                batch_size=BATCH_SIZE,
                convert_to_numpy=True,
                normalize_embeddings=True,
                show_progress_bar=False,
            )
            batch_vecs = batch_vecs.astype(np.float32)

            for j, vec in enumerate(batch_vecs):
                idx = uncached_indices[j]
                text = uncached_texts[j]
                vectors.append((idx, vec))
                self._cache[self._cache_key(text)] = vec

            # 保存增量缓存
            self._save_cache()

        # 按原始顺序重组
        vectors.sort(key=lambda x: x[0])
        return np.stack([v for _, v in vectors])

    # -----------------------------------------------------------------------
    # 相似度计算
    # -----------------------------------------------------------------------
    def compute_scores(self, texts: List[str]) -> np.ndarray:
        """
        计算每条文本与参考概念（默认 "oxidative stress"）的余弦相似度。

        返回
        -------
        np.ndarray, shape (N,)
            相似度分数 ∈ [0, 1]（向量已归一化，等于点积）。
        """
        self._load_model()
        query_vectors = self.encode(texts)
        ref_vector = self._ref_vector
        scores = np.dot(query_vectors, ref_vector)
        return np.clip(scores, 0.0, 1.0)

    # -----------------------------------------------------------------------
    # 分类
    # -----------------------------------------------------------------------
    def classify(
        self,
        texts: List[str],
        threshold: float = DEFAULT_THRESHOLD,
    ) -> List[Tuple[str, float]]:
        """
        对文本列表做语义分类。

        返回
        -------
        List[Tuple[str, float]]
            每个元素为 (category, score)。
            category: "biobert_ox" 或 "biobert_other"。
        """
        scores = self.compute_scores(texts)
        return [
            ("biobert_ox" if s >= threshold else "biobert_other", float(s))
            for s in scores
        ]

    # -----------------------------------------------------------------------
    # 缓存
    # -----------------------------------------------------------------------
    @staticmethod
    def _cache_key(text: str) -> str:
        """基于文本内容的稳定哈希 key。"""
        return hashlib.md5(text.strip().lower().encode()).hexdigest()

    def _load_cache(self) -> None:
        """从磁盘加载缓存。"""
        if self.cache_path.exists():
            try:
                with open(self.cache_path, "rb") as f:
                    data = pickle.load(f)
                if isinstance(data, dict):
                    self._cache = data
            except Exception:
                self._cache = {}

    def _save_cache(self) -> None:
        """保存缓存到磁盘。"""
        self.cache_path.parent.mkdir(parents=True, exist_ok=True)
        with open(self.cache_path, "wb") as f:
            pickle.dump(self._cache, f)

    # -----------------------------------------------------------------------
    # 便捷方法
    # -----------------------------------------------------------------------
    def classify_with_details(
        self,
        texts: List[str],
        threshold: float = DEFAULT_THRESHOLD,
    ) -> Dict[str, List]:
        """
        返回包含详细信息的分类结果（适合直接构造 DataFrame）。

        返回
        -------
        dict with keys: "text", "score", "category", "is_ox"
        """
        scores = self.compute_scores(texts)
        categories = [
            "biobert_ox" if s >= threshold else "biobert_other"
            for s in scores
        ]
        return {
            "text": texts,
            "score": [float(s) for s in scores],
            "category": categories,
            "is_ox": [c == "biobert_ox" for c in categories],
        }


# ---------------------------------------------------------------------------
# 模块级便捷函数（无需实例化）
# ---------------------------------------------------------------------------
_default_classifier: Optional[BioBERTClassifier] = None


def get_classifier() -> BioBERTClassifier:
    """返回全局单例分类器。"""
    global _default_classifier
    if _default_classifier is None:
        _default_classifier = BioBERTClassifier()
    return _default_classifier


def classify_by_biobert(
    texts: List[str],
    threshold: float = DEFAULT_THRESHOLD,
) -> List[Tuple[str, float]]:
    """模块级便捷函数：直接用默认分类器分类。"""
    return get_classifier().classify(texts, threshold=threshold)


def compute_biobert_scores(texts: List[str]) -> np.ndarray:
    """模块级便捷函数：计算语义相似度分数。"""
    return get_classifier().compute_scores(texts)
