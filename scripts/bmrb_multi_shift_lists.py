#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import csv
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from urllib.parse import quote

import requests

API_BASES = ["https://api.bmrb.io", "https://api.bmrb.io/v2"]

SESSION = requests.Session()
SESSION.headers.update({"User-Agent": "bmrb-multi-shift-lists/1.0 (+https://bmrb.io)"})


class ApiError(RuntimeError):
    pass


def _get(url, params=None, timeout=60):
    """GET with small retry/backoff."""
    for attempt in range(4):
        try:
            r = SESSION.get(url, params=params, timeout=timeout)
            if r.status_code == 404:
                # Don’t bother retrying 404s unless we’ll try the other base
                pass
            r.raise_for_status()
            # Some API errors return 200 with HTML text; enforce JSON
            return r.json()
        except requests.HTTPError as e:
            if attempt == 3:
                raise
            time.sleep(0.5 * (2 ** attempt))
        except ValueError:
            # Not JSON
            if attempt == 3:
                raise ApiError(f"Non-JSON response from {url}")
            time.sleep(0.5 * (2 ** attempt))


def fetch_json_across_bases(path, params=None):
    """Try / and /v2 bases for the same path."""
    last_exc = None
    for base in API_BASES:
        url = f"{base}{path}"
        try:
            return _get(url, params=params)
        except Exception as exc:
            last_exc = exc
            continue
    raise ApiError(f"All API base URLs failed for {path}: {last_exc}")


def get_tag_map(tag, database="macromolecules"):
    """
    Returns {bmrb_id (str): [values...]} for the given dictionary tag,
    using /search/get_all_values_for_tag/<tag>.
    """
    # IMPORTANT: URL-encode the tag (leading underscore + dot)
    tag_enc = quote(tag, safe="")
    path = f"/search/get_all_values_for_tag/{tag_enc}"
    return fetch_json_across_bases(path, params={"database": database})


def get_pdb_links_for_entry(bmrb_id):
    """Return [{'pdb_id', 'match_type', 'comment'?}, ...]"""
    path = f"/search/get_pdb_ids_from_bmrb_id/{bmrb_id}"
    try:
        data = fetch_json_across_bases(path)
        return data if isinstance(data, list) else []
    except Exception:
        return []


def main():
    ap = argparse.ArgumentParser(
        description="List BMRB protein entries with ≥2 assigned chemical-shift lists and linked PDB IDs."
    )
    ap.add_argument("-o", "--out", default="bmrb_proteins_multi_shift_lists.csv",
                    help="Output CSV filename")
    ap.add_argument("--all-macromolecules", action="store_true",
                    help="Do NOT filter to proteins; include all macromolecules.")
    ap.add_argument("--max-workers", type=int, default=16,
                    help="Max threads for fetching PDB links (default: 16).")
    args = ap.parse_args()

    # 1) Count assigned shift saveframes per entry
    # Tag per BMRB dictionary: _Assigned_chem_shift_list.Sf_framecode
    shift_sf_map = get_tag_map("_Assigned_chem_shift_list.Sf_framecode")
    shift_sf_list_id = get_tag_map("_Assigned_chem_shift_list.ID")

    # 2) Titles (nice to have)
    title_map = {}
    try:
        title_map = get_tag_map("_Entry.Title")
    except Exception:
        pass

    # 3) Optional: restrict to proteins using _Entity.Polymer_type
    candidate_ids = []
    for bmrb_id, framecodes in (shift_sf_map or {}).items():
        n_lists = len(framecodes) if isinstance(framecodes, list) else (1 if framecodes else 0)
        if n_lists >= 2:
            candidate_ids.append(bmrb_id)

    if not args.all_macromolecules:
        try:
            polymer_type_map = get_tag_map("_Entity.Polymer_type")
        except Exception:
            polymer_type_map = {}
        protein_ids = []
        for bmrb_id in candidate_ids:
            types = polymer_type_map.get(bmrb_id, [])
            if isinstance(types, (list, tuple)):
                is_protein = any(isinstance(t, str) and "polypeptide" in t.lower() for t in types)
            else:
                is_protein = isinstance(types, str) and "polypeptide" in types.lower()
            if is_protein:
                protein_ids.append(bmrb_id)
        candidate_ids = protein_ids

    # 4) Fetch PDB links in parallel
    rows = []
    with ThreadPoolExecutor(max_workers=args.max_workers) as ex:
        future_map = {ex.submit(get_pdb_links_for_entry, bid): bid for bid in candidate_ids}
        for fut in as_completed(future_map):
            bid = future_map[fut]
            pdb_links = fut.result() or []
            pdb_ids = [d.get("pdb_id") for d in pdb_links if isinstance(d, dict) and d.get("pdb_id")]
            match_types = [d.get("match_type") for d in pdb_links if isinstance(d, dict) and d.get("match_type")]
            framecodes = shift_sf_map.get(bid, [])
            n_lists = len(framecodes) if isinstance(framecodes, list) else (1 if framecodes else 0)
            title = ""
            tvals = title_map.get(bid, [])
            if isinstance(tvals, list) and tvals:
                title = tvals[0]
            elif isinstance(tvals, str):
                title = tvals

            rows.append({
                "bmrb_id": bid,
                "title": title,
                "num_shift_lists": n_lists,
                "shift_saveframes": ";".join(framecodes) if isinstance(framecodes, list) else (framecodes or ""),
                "pdb_ids": ";".join(sorted(set(filter(None, pdb_ids)))) if pdb_ids else "",
                "pdb_match_types": ";".join(filter(None, match_types)) if match_types else "",
            })

    # 5) Sort output
    def id_as_int(x):
        try:
            return int(x)
        except Exception:
            return 10**9

    rows.sort(key=lambda r: (-int(r["num_shift_lists"]), id_as_int(r["bmrb_id"])))

    # 6) Write CSV
    fieldnames = ["bmrb_id", "title", "num_shift_lists", "shift_saveframes", "pdb_ids", "pdb_match_types"]
    with open(args.out, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow(r)

    print(f"Wrote {len(rows)} rows to {args.out}", file=sys.stderr)


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"[FATAL] {e}", file=sys.stderr)
        sys.exit(1)
