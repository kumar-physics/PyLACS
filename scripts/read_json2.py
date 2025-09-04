#!/usr/bin/env python3
"""
Fetch all macromolecular entry IDs from BMRB using the BMRB API.

Strategy:
1) Try /v2/list_entries?database=macromolecules (fastest, purpose-built).
2) Fallback: /v2/search/get_all_values_for_tag/_Entry.ID?database=macromolecules
   (returns a dict mapping entry_id -> [values], we take the keys).

Outputs:
- Prints the count and IDs to stdout.
- Optionally writes to a file if --out is provided.

Docs:
- API index: https://api.bmrb.io/v2/
- Databases + endpoints: https://github.com/bmrb-io/BMRB-API
"""

import argparse
import sys
import time
from typing import List, Set
import requests
import os
import json

API_ROOT = "https://api.bmrb.io/v2"
USER_AGENT_HEADER = {"Application": "bmrb-entry-fetcher 1.0"}

def get_with_retry(url: str, params: dict | None = None, max_retries: int = 3, backoff: float = 1.5):
    """GET with simple exponential backoff, returns (ok, response_or_error)."""
    attempt = 0
    while True:
        try:
            r = requests.get(url, params=params, headers=USER_AGENT_HEADER, timeout=30)
            if r.status_code == 200:
                return True, r
            # 403 can be a temporary rate-limit; honor backoff
            if r.status_code in (429, 403, 500, 502, 503, 504):
                attempt += 1
                if attempt > max_retries:
                    return False, f"HTTP {r.status_code}: {r.text[:200]}"
                time.sleep(backoff ** attempt)
            else:
                return False, f"HTTP {r.status_code}: {r.text[:200]}"
        except requests.RequestException as e:
            attempt += 1
            if attempt > max_retries:
                return False, f"Request error: {e}"
            time.sleep(backoff ** attempt)

def fetch_via_list_entries() -> List[str]:
    """
    Use /v2/list_entries?database=macromolecules.
    Expected to return a JSON list of entry IDs (strings).
    """
    url = f"{API_ROOT}/list_entries"
    ok, resp = get_with_retry(url, params={"database": "macromolecules"})
    if not ok:
        raise RuntimeError(f"list_entries failed: {resp}")
    data = resp.json()
    if not isinstance(data, list):
        raise RuntimeError("Unexpected response format from list_entries.")
    # Normalize to strings (some may be ints)
    return sorted(str(x) for x in data)

def fetch_via_get_all_values_for_tag() -> List[str]:
    """
    Fallback using /search/get_all_values_for_tag/_Entry.ID?database=macromolecules
    Returns dict: { "15000": ["15000"], "15001": ["15001"], ... }
    We take the dict keys as the entry IDs.
    """
    url = f"{API_ROOT}/search/get_all_values_for_tag/_Entry.ID"
    ok, resp = get_with_retry(url, params={"database": "macromolecules"})
    if not ok:
        raise RuntimeError(f"get_all_values_for_tag failed: {resp}")
    data = resp.json()
    if not isinstance(data, dict):
        raise RuntimeError("Unexpected response format from get_all_values_for_tag.")
    # Keys are entry IDs; ensure strings
    return sorted(str(k) for k in data.keys())

def fetch_all_macromolecule_ids() -> List[str]:
    """
    Try the primary endpoint first; if it fails, use the fallback.
    """
    try:
        return fetch_via_list_entries()
    except Exception:
        return fetch_via_get_all_values_for_tag()

def main():
    parser = argparse.ArgumentParser(description="Fetch all macromolecular BMRB entry IDs.")
    parser.add_argument("--out", help="Optional path to write IDs (one per line).")
    args = parser.parse_args()

    ids = fetch_all_macromolecule_ids()
    print(f"Total macromolecular BMRB entry IDs: {len(ids)}")
    for eid in ids:
        print(eid)

    if args.out:
        with open(args.out, "w", encoding="utf-8") as f:
            f.write("\n".join(ids))
        print(f"\nWrote IDs to: {args.out}", file=sys.stderr)


def read_json():
    ids = fetch_all_macromolecule_ids()
    for id in ids:
        for method in ["tukey","theilsen","ransac","quantile"]:
            filename = f'/home/nmrbox/kbaskaran/lacs/out/{id}_{method}.json'
            try:
                with open(filename, "r", encoding="utf-8") as f:
                    data = json.load(f)
            except Exception as e:
                print(f"Error reading {filename}: {e}")
            for k in data:
                print (id,k,data[k]['offsets'])
if __name__ == "__main__":
    read_json()
