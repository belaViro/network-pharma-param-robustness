"""Manual DisGeNET v1 endpoint probe.

This script requires DISGENET_API_KEY and performs real network calls only
when executed directly. It is kept outside the default pytest surface.
"""

import os
import time
import requests


BASE = "https://www.disgenet.org/api/v1"


def _headers() -> dict[str, str]:
    api_key = os.environ.get("DISGENET_API_KEY")
    if not api_key:
        raise SystemExit("Set DISGENET_API_KEY before running this manual probe.")
    return {"Authorization": f"Bearer {api_key}", "Accept": "application/json"}


def _print_response(response: requests.Response) -> bool:
    print(f"  status={response.status_code}, type={response.headers.get('Content-Type')}")
    if response.status_code != 200:
        print(f"  -> Text[:200]: {response.text[:200]}")
        return False
    try:
        data = response.json()
        if isinstance(data, list):
            print(f"  -> list len={len(data)}")
            if data:
                print(f"  -> First item: {data[0]}")
                return True
        elif isinstance(data, dict):
            print(f"  -> JSON keys: {list(data.keys())}")
            print(f"  -> Data preview: {str(data)[:500]}")
            return bool(data)
    except ValueError:
        print(f"  -> Text[:300]: {response.text[:300]}")
    return False


def main() -> None:
    headers = _headers()

    print("=== Call 1: entity/disease ===")
    for param_name in ["disease_name", "disease", "query", "name"]:
        url = f"{BASE}/entity/disease?{param_name}=Oxidative%20Stress"
        print(f"  {param_name}:")
        if _print_response(requests.get(url, headers=headers, timeout=15)):
            break
        time.sleep(1)

    print("\n=== Call 2: gda/summary (by disease name) ===")
    for param_name in ["disease_name", "disease", "diseaseId", "disease_id"]:
        url = f"{BASE}/gda/summary?{param_name}=Oxidative%20Stress"
        print(f"  {param_name}:")
        if _print_response(requests.get(url, headers=headers, timeout=15)):
            break
        time.sleep(1)

    print("\n=== Call 3: gda/summary (by CUI) ===")
    for disease_id in ["C0242606", "C0029682"]:
        url = f"{BASE}/gda/summary?diseaseId={disease_id}"
        print(f"  {disease_id}:")
        if _print_response(requests.get(url, headers=headers, timeout=15)):
            break
        time.sleep(1)


if __name__ == "__main__":
    main()
