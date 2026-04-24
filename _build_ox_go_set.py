"""
用 goatools 构建氧化应激相关 GO term 集合
以 GO:0006979 (response to oxidative stress) 为根，获取所有后代节点
"""
import os
import requests
from goatools.obo_parser import GODag

OBO_URL = "http://purl.obolibrary.org/obo/go/go-basic.obo"
OBO_PATH = "src/data/go-basic.obo"

def download_obo():
    if os.path.exists(OBO_PATH):
        print(f"OBO file already exists: {OBO_PATH}")
        return
    os.makedirs(os.path.dirname(OBO_PATH), exist_ok=True)
    print(f"Downloading GO OBO file from {OBO_URL}...")
    headers = {"User-Agent": "Mozilla/5.0 (research script)"}
    r = requests.get(OBO_URL, headers=headers, timeout=120)
    r.raise_for_status()
    with open(OBO_PATH, 'wb') as f:
        f.write(r.content)
    print(f"Saved to {OBO_PATH} ({len(r.content)} bytes)")

def build_ox_stress_set():
    download_obo()
    print("Loading GO DAG (this may take a few seconds)...")
    dag = GODag(OBO_PATH)
    
    root_id = "GO:0006979"
    if root_id not in dag:
        print(f"ERROR: {root_id} not found in GO DAG!")
        return set()
    
    root = dag[root_id]
    print(f"\nRoot term: {root_id} - {root.name}")
    
    # 获取所有后代（descendants）
    descendants = root.get_all_children()
    descendants.add(root_id)
    
    print(f"Total oxidative-stress-related GO terms: {len(descendants)}")
    
    # 展示一些示例
    print("\nSample terms (first 30):")
    for go_id in sorted(list(descendants))[:30]:
        term = dag.get(go_id)
        name = term.name if term else "N/A"
        print(f"  {go_id}: {name}")
    
    return descendants, dag

if __name__ == "__main__":
    ox_set, dag = build_ox_stress_set()
    
    # 保存到文件供后续使用
    out_path = "src/data/oxidative_stress_go_terms.txt"
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, 'w', encoding='utf-8') as f:
        for go_id in sorted(ox_set):
            term = dag.get(go_id)
            name = term.name if term else "N/A"
            f.write(f"{go_id}\t{name}\n")
    print(f"\nSaved to {out_path}")
