"""Manual DisGeNET endpoint probe.

This file is intentionally not a pytest test. It performs real network calls
only when executed directly and reads the API token from DISGENET_API_KEY.
"""

import os
import requests


URLS = [
    "https://www.disgenet.com/api/disease/search/Oxidative%20Stress",
    "https://www.disgenet.org/api/disease/search/Oxidative%20Stress",
    "https://api.disgenet.com/disease/search/Oxidative%20Stress",
]


def main() -> None:
    api_key = os.environ.get("DISGENET_API_KEY")
    if not api_key:
        raise SystemExit("Set DISGENET_API_KEY before running this manual probe.")

    headers = {"Authorization": f"Bearer {api_key}", "Accept": "application/json"}
    for url in URLS:
        try:
            response = requests.get(url, headers=headers, timeout=10)
            print(f"URL: {url}")
            print(f"  Status: {response.status_code}")
            print(f"  Content-Type: {response.headers.get('Content-Type')}")
            print(f"  Text[:300]: {response.text[:300]}")
            print()
        except Exception as exc:
            print(f"URL: {url}")
            print(f"  ERROR: {exc}")
            print()


if __name__ == "__main__":
    main()
