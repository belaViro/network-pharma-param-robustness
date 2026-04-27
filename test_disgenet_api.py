import requests

# 尝试找到正确的 API 端点
urls = [
    "https://www.disgenet.org/swagger.json",
    "https://www.disgenet.org/api/swagger.json",
    "https://www.disgenet.org/api/v1/",
    "https://www.disgenet.org/api/v1/disease/search/Oxidative%20Stress",
]

for url in urls:
    try:
        r = requests.get(url, timeout=10)
        ct = r.headers.get("Content-Type", "unknown")
        print(f"{url}: status={r.status_code}, type={ct}")
    except Exception as e:
        print(f"{url}: ERROR {e}")
