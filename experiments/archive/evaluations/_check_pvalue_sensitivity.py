"""
验证核心组合在不同P-value标准下的富集结果差异
核心组合: P=0.05, STRING=0.4, TopK=30
P-value标准: 0.001, 0.01, 0.05
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))

import pandas as pd
import networkx as nx
from src.config import SWISS_TARGETS_DIR, BACKGROUND_GENES
from src.enrichment.enrichment_analysis import run_enrichment_analysis

# 1. 加载数据
all_predictions = []
for csv_file in SWISS_TARGETS_DIR.glob("*.csv"):
    df = pd.read_csv(csv_file)
    if "Target" in df.columns and "Common name" in df.columns and "Probability*" in df.columns:
        all_predictions.append(df[["Common name", "Probability*"]].rename(
            columns={"Common name": "Gene", "Probability*": "Probability"}
        ))

pred_df = pd.concat(all_predictions, ignore_index=True).drop_duplicates(subset=["Gene"])

# 2. 筛选靶点 (Prob>=0.05)
filtered_targets = pred_df[pred_df["Probability"] >= 0.05]["Gene"].tolist()
print(f"Prob>=0.05 靶点数: {len(filtered_targets)}")

# 3. 构建PPI网络 (STRING>=0.4)
string_edges = pd.read_csv("src/outputs/shared/ppi_network/string_interactions_interactions.csv")
filtered_edges = string_edges[string_edges["combined_score"] >= 0.4]

G = nx.Graph()
for _, row in filtered_edges.iterrows():
    G.add_edge(row["protein1"], row["protein2"], weight=row["combined_score"])

valid_nodes = [n for n in G.nodes() if n in filtered_targets]
G = G.subgraph(valid_nodes).copy()
print(f"PPI网络: {G.number_of_nodes()} 节点, {G.number_of_edges()} 边")

# 4. 取Top-30 Hub
degree_sorted = sorted(G.degree(), key=lambda x: x[1], reverse=True)
topk = [g for g, d in degree_sorted[:30]]
print(f"\nTop-30 Hub靶点: {', '.join(topk)}")

# 5. 三个P-value标准下的富集分析
PVALUES = [0.001, 0.01, 0.05]

core_pathways = {
    "MAPK signaling pathway": "KEGG",
    "PI3K-Akt signaling pathway": "KEGG",
    "Rap1 signaling pathway": "KEGG",
    "Ras signaling pathway": "KEGG",
    "AGE-RAGE signaling pathway in diabetic complications": "KEGG",
    "IL-17 signaling pathway": "KEGG",
    "Fluid shear stress and atherosclerosis": "KEGG",
    "Pathways in cancer": "KEGG",
    "Proteoglycans in cancer": "KEGG",
}

results = []

for pval in PVALUES:
    print(f"\n{'='*60}")
    print(f"P-value cutoff = {pval}")
    print(f"{'='*60}")
    
    enrichment = run_enrichment_analysis(
        gene_list=topk,
        pvalue_cutoff=pval,
        top_n=20,
        background=BACKGROUND_GENES
    )
    
    total_sig = 0
    if enrichment["status"] == "成功" and "significant" in enrichment:
        for source, df_sig in enrichment["significant"].items():
            if isinstance(df_sig, pd.DataFrame) and len(df_sig) > 0:
                total_sig += len(df_sig)
    
    print(f"总显著通路数: {total_sig}")
    print(f"\n核心通路在该标准下的显著性:")
    
    for pathway_name, expected_source in core_pathways.items():
        found = False
        pvalue = None
        if enrichment["status"] == "成功" and "significant" in enrichment:
            for source, df_sig in enrichment["significant"].items():
                if isinstance(df_sig, pd.DataFrame) and len(df_sig) > 0:
                    for _, row in df_sig.iterrows():
                        if row.get("Term", "") == pathway_name:
                            found = True
                            pvalue = row.get("Adjusted P-value", row.get("P-value", None))
                            break
        
        status = "[SIG]" if found else "[NS] "
        pv_str = f" (P={pvalue:.2e})" if pvalue else ""
        print(f"  {pathway_name:<55} {status}{pv_str}")
        
        results.append({
            "pvalue_cutoff": pval,
            "pathway": pathway_name,
            "significant": found,
            "actual_pvalue": pvalue,
        })

# 6. 汇总对比
print(f"\n{'='*60}")
print("汇总对比表")
print(f"{'='*60}")

results_df = pd.DataFrame(results)
pivot = results_df.pivot(index="pathway", columns="pvalue_cutoff", values="significant")
print(pivot.to_string())

print(f"\n结论:")
core_always_sig = all(pivot.loc["MAPK signaling pathway"]) and \
                  all(pivot.loc["PI3K-Akt signaling pathway"])
if core_always_sig:
    print("  核心通路（MAPK/PI3K-Akt）在所有三个P-value标准下均显著")
    print("  → P_pvalue扫描的实际信息增益有限")
else:
    print("  部分核心通路在严格标准下消失 → P_pvalue扫描有价值")
