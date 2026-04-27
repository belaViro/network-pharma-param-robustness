import requests

api_key = "62d0ad0a-da18-4761-8f59-77b2fc207661"
headers = {"Authorization": f"Bearer {api_key}", "Accept": "application/json"}

urls = [
    "https://www.disgenet.com/api/disease/search/Oxidative%20Stress",
    "https://www.disgenet.org/api/disease/search/Oxidative%20Stress",
    "https://api.disgenet.com/disease/search/Oxidative%20Stress",
]

for url in urls:
    try:
        r = requests.get(url, headers=headers, timeout=10)
        print(f"URL: {url}")
        print(f"  Status: {r.status_code}")
        print(f"  Content-Type: {r.headers.get('Content-Type')}")
        print(f"  Text[:300]: {r.text[:300]}")
        print()
    except Exception as e:
        print(f"URL: {url}")
        print(f"  ERROR: {e}")
        print()
