"""
每个参数组合单独富集 vs Q1稳健靶点富集 横向对比
目的：量化单组合结果的波动，证明稳健靶点的价值
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))

import pandas as pd
import numpy as np
import networkx as nx
from itertools import product
import json
import hashlib
import re

from src.config import SWISS_TARGETS_DIR, BACKGROUND_GENES
from src.enrichment.enrichment_analysis import run_enrichment_analysis

# ============================================================================
# 0. 缓存
# ============================================================================
CACHE_DIR = Path("src/outputs/sensitivity/enrichment_cache")
CACHE_DIR.mkdir(parents=True, exist_ok=True)

def get_cache_key(gene_list, pvalue_cutoff):
    genes_str = ";".join(sorted(gene_list))
    return hashlib.md5(f"{genes_str}|{pvalue_cutoff}".encode()).hexdigest()[:12]

def load_cached(cache_key):
    f = CACHE_DIR / f"{cache_key}.json"
    if f.exists():
        with open(f, "r") as fh:
            data = json.load(fh)
            return {k: pd.DataFrame(v) if isinstance(v, list) else v for k, v in data.items()}
    return None

def save_cached(cache_key, result):
    f = CACHE_DIR / f"{cache_key}.json"
    serializable = {k: (v.to_dict(orient="records") if isinstance(v, pd.DataFrame) else v) for k, v in result.items()}
    with open(f, "w") as fh:
        json.dump(serializable, fh)

def cached_enrichment(gene_list, pvalue_cutoff, background=20500):
    ck = get_cache_key(gene_list, pvalue_cutoff)
    cached = load_cached(ck)
    if cached is not None:
        return {"status": "成功", "significant": cached}
    enrichment = run_enrichment_analysis(gene_list=gene_list, pvalue_cutoff=pvalue_cutoff, top_n=20, background=background)
    if enrichment["status"] == "成功" and "significant" in enrichment:
        save_cached(ck, enrichment["significant"])
    return enrichment

# ============================================================================
# 1. GO集合 + 分类
# ============================================================================
GO_OX_FILE = Path("src/data/oxidative_stress_go_terms.txt")
GO_OX_IDS = set()
if GO_OX_FILE.exists():
    with open(GO_OX_FILE, 'r', encoding='utf-8') as f:
        for line in f:
            go_id = line.strip().split('\t')[0]
            GO_OX_IDS.add(go_id)

CROSS_TALK_PATHWAYS = {
    "AGE-RAGE signaling pathway in diabetic complications",
    "IL-17 signaling pathway",
    "Ferroptosis",
    "Glutathione metabolism",
    "Peroxisome",
    "Fluid shear stress and atherosclerosis",
    "NOD-like receptor signaling pathway",
    "TNF signaling pathway",
    "Toll-like receptor signaling pathway",
    "NF-kappa B signaling pathway",
}

BACKBONE_KEYWORDS = [
    'MAPK', 'PI3K-Akt', 'Rap1', 'Ras signaling',
    'receptor tyrosine kinase', 'cytokine-mediated',
    'calcium signaling', 'cAMP signaling',
    'chemokine signaling', 'JAK-STAT',
    'EGFR tyrosine kinase inhibitor resistance',
]

def classify_pathway(term: str) -> str:
    term_lower = term.lower()
    match = re.search(r'\(GO:(\d+)\)', term)
    if match:
        go_id = f"GO:{match.group(1)}"
        if go_id in GO_OX_IDS:
            return "direct_ox"
    if term in CROSS_TALK_PATHWAYS:
        return "cross_talk"
    for kw in BACKBONE_KEYWORDS:
        if re.search(r'\b' + re.escape(kw.lower()) + r'\b', term_lower):
            return "signaling_backbone"
    return "other"

# ============================================================================
# 2. 数据加载
# ============================================================================
PROB_THRESHOLDS = [0.0, 0.05, 0.1]
STRING_THRESHOLDS = [0.3, 0.4, 0.5, 0.7]
PVALUE_CUTOFFS = [0.01]
TOP_K_VALUES = [25, 30, 40]
PARAM_COMBINATIONS = list(product(PROB_THRESHOLDS, STRING_THRESHOLDS, PVALUE_CUTOFFS, TOP_K_VALUES))

all_predictions = []
for csv_file in SWISS_TARGETS_DIR.glob("*.csv"):
    df = pd.read_csv(csv_file)
    if "Target" in df.columns and "Common name" in df.columns and "Probability*" in df.columns:
        all_predictions.append(df[["Common name", "Probability*"]].rename(
            columns={"Common name": "Gene", "Probability*": "Probability"}
        ))
pred_df = pd.concat(all_predictions, ignore_index=True).drop_duplicates(subset=["Gene"])
string_edges = pd.read_csv("src/outputs/shared/ppi_network/string_interactions_interactions.csv")

# ============================================================================
# 3. Q1稳健靶点（从mvp_v2结果读取）
# ============================================================================
q1_df = pd.read_csv("src/outputs/sensitivity/mvp_v2_target_robustness.csv")
q1_targets = q1_df[q1_df["quadrant"].str.contains("Q1", na=False)]["gene"].tolist()
print(f"Q1稳健靶点: {len(q1_targets)} 个")
print(f"  {', '.join(q1_targets)}")

# ============================================================================
# 4. 对Q1稳健靶点做富集
# ============================================================================
q1_enrichment = cached_enrichment(q1_targets, 0.01, BACKGROUND_GENES)
q1_pathways = {}  # pathway -> {logp, ef, category}
if q1_enrichment["status"] == "成功" and "significant" in q1_enrichment:
    for source, df_sig in q1_enrichment["significant"].items():
        if isinstance(df_sig, pd.DataFrame) and len(df_sig) > 0:
            for _, row in df_sig.iterrows():
                term = str(row.get("Term", ""))
                pval = row.get("Adjusted P-value", 1.0)
                logp = -np.log10(max(pval, 1e-300))
                ef = row.get("EF", 0)
                cat = classify_pathway(term)
                q1_pathways[term] = {"logp": logp, "ef": ef, "category": cat, "source": source}

print(f"\nQ1稳健靶点显著通路: {len(q1_pathways)} 条")
q1_backbone = [p for p, info in q1_pathways.items() if info["category"] == "signaling_backbone"]
q1_crosstalk = [p for p, info in q1_pathways.items() if info["category"] == "cross_talk"]
print(f"  主干通路: {len(q1_backbone)} 条")
print(f"  交叉调控: {len(q1_crosstalk)} 条")

# ============================================================================
# 5. 对每个参数组合跑富集
# ============================================================================
print("\n" + "=" * 70)
print("每个参数组合单独富集分析")
print("=" * 70)

combo_results = []
combo_pathway_records = []

for idx, (prob_th, string_th, pval_th, top_k) in enumerate(PARAM_COMBINATIONS):
    combo_str = f"P{prob_th}_S{string_th}_V{pval_th}_K{top_k}"
    filtered_targets = pred_df[pred_df["Probability"] >= prob_th]["Gene"].tolist()
    
    if len(filtered_targets) < 10:
        continue
    
    filtered_edges = string_edges[string_edges["combined_score"] >= string_th]
    if len(filtered_edges) == 0:
        continue
    
    G = nx.Graph()
    for _, row in filtered_edges.iterrows():
        G.add_edge(row["protein1"], row["protein2"], weight=row["combined_score"])
    
    valid_nodes = [n for n in G.nodes() if n in filtered_targets]
    G = G.subgraph(valid_nodes).copy()
    if G.number_of_nodes() < 10:
        continue
    
    degree_sorted = sorted(G.degree(), key=lambda x: x[1], reverse=True)
    topk_genes = [g for g, d in degree_sorted[:top_k]]
    if len(topk_genes) < 5:
        continue
    
    # 富集
    enrichment = cached_enrichment(topk_genes, pval_th, BACKGROUND_GENES)
    
    sig_pathways = {}
    total_sig = 0
    avg_logp = 0.0
    max_logp = 0.0
    backbone_count = 0
    crosstalk_count = 0
    
    if enrichment["status"] == "成功" and "significant" in enrichment:
        for source, df_sig in enrichment["significant"].items():
            if isinstance(df_sig, pd.DataFrame) and len(df_sig) > 0:
                total_sig += len(df_sig)
                for _, row in df_sig.iterrows():
                    term = str(row.get("Term", ""))
                    pval = row.get("Adjusted P-value", 1.0)
                    logp = -np.log10(max(pval, 1e-300))
                    ef = row.get("EF", 0)
                    cat = classify_pathway(term)
                    
                    sig_pathways[term] = {"logp": logp, "ef": ef, "category": cat, "source": source}
                    avg_logp += logp
                    if logp > max_logp:
                        max_logp = logp
                    if cat == "signaling_backbone":
                        backbone_count += 1
                    elif cat == "cross_talk":
                        crosstalk_count += 1
                    
                    combo_pathway_records.append({
                        "combo": combo_str,
                        "prob_th": prob_th,
                        "string_th": string_th,
                        "top_k": top_k,
                        "pathway": term,
                        "logp": logp,
                        "ef": ef,
                        "category": cat,
                        "source": source,
                    })
        
        if total_sig > 0:
            avg_logp /= total_sig
    
    # 与Q1的重叠度
    q1_overlap = len(set(sig_pathways.keys()) & set(q1_pathways.keys()))
    q1_overlap_pct = q1_overlap / len(q1_pathways) * 100 if q1_pathways else 0
    
    combo_results.append({
        "combo": combo_str,
        "prob_th": prob_th,
        "string_th": string_th,
        "top_k": top_k,
        "n_targets": len(filtered_targets),
        "topk_size": len(topk_genes),
        "n_sig": total_sig,
        "avg_logp": avg_logp,
        "max_logp": max_logp,
        "backbone_count": backbone_count,
        "crosstalk_count": crosstalk_count,
        "q1_overlap": q1_overlap,
        "q1_overlap_pct": q1_overlap_pct,
    })
    
    print(f"  [{idx+1:2d}/36] {combo_str}: 靶点{len(filtered_targets):3d} TopK{len(topk_genes):2d} "
          f"通路{total_sig:3d} 主干{backbone_count:2d} 交叉{crosstalk_count:2d} "
          f"avg_logp={avg_logp:.1f} Q1重叠={q1_overlap_pct:.0f}%")

combo_df = pd.DataFrame(combo_results)

# ============================================================================
# 6. 波动分析
# ============================================================================
print("\n" + "=" * 70)
print("单组合结果波动分析")
print("=" * 70)

print(f"\n显著通路数: 均值={combo_df['n_sig'].mean():.1f}, 标准差={combo_df['n_sig'].std():.1f}, "
      f"范围=[{combo_df['n_sig'].min()}, {combo_df['n_sig'].max()}]")
print(f"平均-logP:  均值={combo_df['avg_logp'].mean():.1f}, 标准差={combo_df['avg_logp'].std():.1f}, "
      f"范围=[{combo_df['avg_logp'].min():.1f}, {combo_df['avg_logp'].max():.1f}]")
print(f"最大-logP:  均值={combo_df['max_logp'].mean():.1f}, 标准差={combo_df['max_logp'].std():.1f}, "
      f"范围=[{combo_df['max_logp'].min():.1f}, {combo_df['max_logp'].max():.1f}]")
print(f"主干通路数: 均值={combo_df['backbone_count'].mean():.1f}, 标准差={combo_df['backbone_count'].std():.1f}, "
      f"范围=[{combo_df['backbone_count'].min()}, {combo_df['backbone_count'].max()}]")
print(f"Q1重叠度:   均值={combo_df['q1_overlap_pct'].mean():.1f}%, 标准差={combo_df['q1_overlap_pct'].std():.1f}%, "
      f"范围=[{combo_df['q1_overlap_pct'].min():.0f}%, {combo_df['q1_overlap_pct'].max():.0f}%]")

# ============================================================================
# 7. 通路级波动：哪些通路"时隐时现"
# ============================================================================
print("\n" + "=" * 70)
print("通路级稳定性：哪些通路在36个组合中时隐时现？")
print("=" * 70)

cpr = pd.DataFrame(combo_pathway_records)
pathway_stability = []
for pathway, group in cpr.groupby("pathway"):
    freq = len(group) / len(combo_results)
    avg_logp = group["logp"].mean()
    logp_std = group["logp"].std()
    category = group["category"].iloc[0]
    
    # 是否在Q1中出现
    in_q1 = pathway in q1_pathways
    q1_logp = q1_pathways[pathway]["logp"] if in_q1 else None
    
    pathway_stability.append({
        "pathway": pathway,
        "category": category,
        "frequency": freq,
        "count": len(group),
        "avg_logp": avg_logp,
        "logp_std": logp_std,
        "in_q1": in_q1,
        "q1_logp": q1_logp,
    })

pathway_stab_df = pd.DataFrame(pathway_stability).sort_values("frequency", ascending=False)

print(f"\n{'通路':<55} {'类别':<10} {'出现频率':<8} {'平均-logP':<10} {'Q1中有?':<8}")
print("-" * 100)
for _, row in pathway_stab_df.head(20).iterrows():
    cat_short = {"direct_ox": "OX", "cross_talk": "CT", "signaling_backbone": "BB", "other": "OT"}.get(row["category"], row["category"])
    q1_mark = "YES" if row["in_q1"] else "no"
    print(f"{row['pathway']:<53} {cat_short:<10} {row['frequency']:<8.2f} {row['avg_logp']:<10.1f} {q1_mark:<8}")

# ============================================================================
# 8. 核心问题：Q1稳健靶点 vs 单组合平均
# ============================================================================
print("\n" + "=" * 70)
print("核心对比：Q1稳健靶点 vs 36个单组合的平均表现")
print("=" * 70)

q1_avg_logp = np.mean([info["logp"] for info in q1_pathways.values()])
q1_max_logp = max([info["logp"] for info in q1_pathways.values()])
q1_n_sig = len(q1_pathways)
q1_n_bb = len([p for p, info in q1_pathways.items() if info["category"] == "signaling_backbone"])
q1_n_ct = len([p for p, info in q1_pathways.items() if info["category"] == "cross_talk"])

print(f"\n{'指标':<25} {'Q1稳健靶点':<15} {'单组合均值':<15} {'单组合范围':<20}")
print("-" * 75)
print(f"{'显著通路数':<25} {q1_n_sig:<15} {combo_df['n_sig'].mean():<15.1f} [{combo_df['n_sig'].min()}, {combo_df['n_sig'].max()}]")
print(f"{'平均-logP':<25} {q1_avg_logp:<15.1f} {combo_df['avg_logp'].mean():<15.1f} [{combo_df['avg_logp'].min():.1f}, {combo_df['avg_logp'].max():.1f}]")
print(f"{'最大-logP':<25} {q1_max_logp:<15.1f} {combo_df['max_logp'].mean():<15.1f} [{combo_df['max_logp'].min():.1f}, {combo_df['max_logp'].max():.1f}]")
print(f"{'主干通路数':<25} {q1_n_bb:<15} {combo_df['backbone_count'].mean():<15.1f} [{combo_df['backbone_count'].min()}, {combo_df['backbone_count'].max()}]")
print(f"{'交叉调控通路数':<25} {q1_n_ct:<15} {combo_df['crosstalk_count'].mean():<15.1f} [{combo_df['crosstalk_count'].min()}, {combo_df['crosstalk_count'].max()}]")

# ============================================================================
# 9. 波动最大的通路（不稳定的发现）
# ============================================================================
print("\n" + "=" * 70)
print("波动最大的通路（'时隐时现'的发现）")
print("=" * 70)

unstable = pathway_stab_df[
    (pathway_stab_df["frequency"] >= 0.3) & 
    (pathway_stab_df["frequency"] <= 0.8)
].sort_values("logp_std", ascending=False)

print(f"\n出现频率30%-80%的通路（不稳定区域）共 {len(unstable)} 条:")
print(f"{'通路':<55} {'频率':<8} {'-logP波动':<10} {'Q1中有?'}")
print("-" * 85)
for _, row in unstable.head(15).iterrows():
    q1_mark = "YES" if row["in_q1"] else "no"
    print(f"{row['pathway']:<53} {row['frequency']:<8.2f} {row['logp_std']:<10.1f} {q1_mark}")

# ============================================================================
# 10. 保存结果
# ============================================================================
combo_df.to_csv("src/outputs/sensitivity/per_combo_summary.csv", index=False)
cpr.to_csv("src/outputs/sensitivity/per_combo_pathways.csv", index=False)
pathway_stab_df.to_csv("src/outputs/sensitivity/pathway_stability.csv", index=False)

print(f"\n结果已保存:")
print(f"  - per_combo_summary.csv: 36个组合汇总")
print(f"  - per_combo_pathways.csv: 每个组合的通路详情")
print(f"  - pathway_stability.csv: 通路级稳定性")
