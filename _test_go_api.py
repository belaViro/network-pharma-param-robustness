import requests
import re

def fetch_ox_stress_descendants():
    """从 GO API 获取 GO:0006979 (response to oxidative stress) 的所有后代"""
    url = "https://api.geneontology.org/api/ontology/term/GO:0006979/graph"
    params = {"relationship_type": "is_a", "descendants": "true"}
    
    print("Querying GO API...")
    r = requests.get(url, params=params, timeout=60)
    print(f"Status: {r.status_code}")
    
    data = r.json()
    nodes = data.get("nodes", [])
    edges = data.get("edges", [])
    
    print(f"Total nodes in graph: {len(nodes)}")
    print(f"Total edges in graph: {len(edges)}")
    
    # 提取所有 GO IDs
    go_ids = {node["id"] for node in nodes}
    go_ids.add("GO:0006979")
    
    # 展示一些示例
    print(f"\n氧化应激相关 GO terms 总数: {len(go_ids)}")
    print("示例（前20条）:")
    for node in nodes[:20]:
        print(f"  {node['id']}: {node.get('label', 'N/A')}")
    
    return go_ids

def extract_go_id(term: str) -> str | None:
    match = re.search(r'\(GO:(\d+)\)', term)
    if match:
        return f"GO:{match.group(1)}"
    return None

# 测试
if __name__ == "__main__":
    ox_ids = fetch_ox_stress_descendants()
    
    # 测试咱们实际数据中的一些 GO terms
    test_terms = [
        "positive regulation of phosphorylation (GO:0042327)",
        "transmembrane receptor protein tyrosine kinase signaling pathway (GO:0007169)",
        "cellular response to oxidative stress (GO:0034599)",
        "response to oxidative stress (GO:0006979)",
        "positive regulation of MAP kinase activity (GO:0043406)",
        "protein tyrosine kinase activity (GO:0004713)",
        "membrane raft (GO:0045121)",
    ]
    
    print("\n测试判断:")
    for term in test_terms:
        go_id = extract_go_id(term)
        is_ox = go_id in ox_ids if go_id else False
        status = "OX" if is_ox else "NOT OX"
        print(f"  [{status}] {term}")
