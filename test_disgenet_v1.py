import requests, json, time

API_KEY = "62d0ad0a-da18-4761-8f59-77b2fc207661"
HEADERS = {
    "Authorization": f"Bearer {API_KEY}",
    "Accept": "application/json",
}
BASE = "https://www.disgenet.org/api/v1"

# Call 1: 搜索疾病 ID（尝试多种参数名）
print("=== Call 1: entity/disease ===")
for param_name in ["disease_name", "disease", "query", "name"]:
    url = f"{BASE}/entity/disease?{param_name}=Oxidative%20Stress"
    r = requests.get(url, headers=HEADERS, timeout=15)
    print(f"  {param_name}: status={r.status_code}, type={r.headers.get('Content-Type')}")
    if r.status_code == 200:
        try:
            data = r.json()
            print(f"  -> JSON keys: {list(data.keys()) if isinstance(data, dict) else 'list len=' + str(len(data))}")
            if isinstance(data, list) and len(data) > 0:
                print(f"  -> First item: {data[0]}")
                break
            elif isinstance(data, dict) and data:
                print(f"  -> Data preview: {str(data)[:500]}")
                break
        except:
            print(f"  -> Text[:300]: {r.text[:300]}")
    else:
        print(f"  -> Text[:200]: {r.text[:200]}")
    time.sleep(1)

# Call 2: 直接尝试 gda/summary 用疾病名称
print("\n=== Call 2: gda/summary (by disease name) ===")
for param_name in ["disease_name", "disease", "diseaseId", "disease_id"]:
    url = f"{BASE}/gda/summary?{param_name}=Oxidative%20Stress"
    r = requests.get(url, headers=HEADERS, timeout=15)
    print(f"  {param_name}: status={r.status_code}")
    if r.status_code == 200:
        try:
            data = r.json()
            print(f"  -> JSON keys: {list(data.keys()) if isinstance(data, dict) else 'list len=' + str(len(data))}")
            if isinstance(data, list) and len(data) > 0:
                print(f"  -> First item keys: {list(data[0].keys())}")
                print(f"  -> Sample: {data[0]}")
                break
        except:
            print(f"  -> Text[:300]: {r.text[:300]}")
    else:
        print(f"  -> Text[:200]: {r.text[:200]}")
    time.sleep(1)

# Call 3: 尝试 disease ID (CUI)
print("\n=== Call 3: gda/summary (by CUI) ===")
for disease_id in ["C0242606", "C0029682"]:
    url = f"{BASE}/gda/summary?diseaseId={disease_id}"
    r = requests.get(url, headers=HEADERS, timeout=15)
    print(f"  {disease_id}: status={r.status_code}")
    if r.status_code == 200:
        try:
            data = r.json()
            print(f"  -> JSON keys: {list(data.keys()) if isinstance(data, dict) else 'list len=' + str(len(data))}")
            if isinstance(data, list) and len(data) > 0:
                print(f"  -> First item keys: {list(data[0].keys())}")
                print(f"  -> Total associations: {len(data)}")
                break
        except:
            print(f"  -> Text[:300]: {r.text[:300]}")
    else:
        print(f"  -> Text[:200]: {r.text[:200]}")
    time.sleep(1)
