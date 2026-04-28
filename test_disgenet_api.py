"""Manual DisGeNET API discovery probe.

This module has no import-time side effects so pytest collection stays offline.
Run it directly only when network access is intended.
"""

import requests


URLS = [
    "https://www.disgenet.org/swagger.json",
    "https://www.disgenet.org/api/swagger.json",
    "https://www.disgenet.org/api/v1/",
    "https://www.disgenet.org/api/v1/disease/search/Oxidative%20Stress",
]


def main() -> None:
    for url in URLS:
        try:
            response = requests.get(url, timeout=10)
            content_type = response.headers.get("Content-Type", "unknown")
            print(f"{url}: status={response.status_code}, type={content_type}")
        except Exception as exc:
            print(f"{url}: ERROR {exc}")


if __name__ == "__main__":
    main()
