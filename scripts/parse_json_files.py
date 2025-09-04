#!/usr/bin/env python3
"""
Parse .json files inside folders of a directory.

Usage:
  python parse_json_in_folders.py /path/to/dir
  python parse_json_in_folders.py /path/to/dir --recursive
  python parse_json_in_folders.py /path/to/dir --include-root
"""

from __future__ import annotations
import argparse
import json
from pathlib import Path
from typing import Any, Dict, Iterable, Tuple
import numpy as np
import plotly.express as px
import pandas as pd
def iter_json_files(
    root: Path,
    recursive: bool = False,
    include_root: bool = False,
) -> Iterable[Path]:
    """
    Yield .json files under each folder in `root`.

    - If include_root is False (default), only files inside subdirectories are yielded.
    - If recursive is True, descends into nested subfolders as well.
    """
    if not root.is_dir():
        raise NotADirectoryError(f"{root} is not a directory")

    if include_root:
        # Optionally include JSONs directly in the root dir
        for p in root.glob("*.json"):
            if p.is_file():
                yield p

    if recursive:
        for p in root.rglob("*.json"):
            # Skip those already yielded from root if include_root, otherwise rglob covers everything
            if include_root and p.parent == root:
                continue
            if p.is_file():
                yield p
    else:
        # One level deep: iterate immediate subdirectories only
        for sub in (p for p in root.iterdir() if p.is_dir()):
            for jf in sub.glob("*.json"):
                if jf.is_file():
                    yield jf


def parse_json_file(path: Path) -> Tuple[Path, Dict[str, Any] | None, str | None]:
    """
    Try to load a JSON file and return (path, data, error_message).
    On success, error_message is None. On failure, data is None.
    """
    try:
        with path.open("r", encoding="utf-8") as f:
            data = json.load(f)
        return path, data, None
    except json.JSONDecodeError as e:
        return path, None, f"JSON decode error at line {e.lineno} col {e.colno}: {e.msg}"
    except OSError as e:
        return path, None, f"OS error: {e.strerror} ({e.filename})"
    except Exception as e:
        return path, None, f"Unexpected error: {type(e).__name__}: {e}"


def process_json(data: Dict[str, Any], path: Path) -> None:
    """
    Replace with your custom logic.
    For demo: print file path and top-level keys (up to 10).
    """


    keys = list(data.keys())
    preview = ", ".join(keys[:10]) + (" ..." if len(keys) > 10 else "")
    if len(preview):
        #print(f"OK  | {path} | keys: [{preview}]")
        entry_id = path.stem.split("_")[0]
        method = path.stem.split("_")[-1]
        for list_id in data:
            for atom in data[list_id]['offsets']:
                print (f'{entry_id},{method},{atom},{data[list_id]['offsets'][atom]}')



def main() -> None:
    ap = argparse.ArgumentParser(description="Parse JSON files in subfolders.")
    ap.add_argument("directory", type=Path, help="Root directory to scan")
    ap.add_argument(
        "--recursive",
        action="store_true",
        help="Recurse into nested subfolders (default: only immediate subfolders)",
    )
    ap.add_argument(
        "--include-root",
        action="store_true",
        help="Also parse JSON files directly in the root directory",
    )
    args = ap.parse_args()

    total = 0
    ok = 0
    errors = 0
    dict_data = {'tukey': {'ca': [], 'cb': [], 'c': [], 'n': [], 'ha': []},
                 'ransac': {'ca': [], 'cb': [], 'c': [], 'n': [], 'ha': []},
                 'quantile': {'ca': [], 'cb': [], 'c': [], 'n': [], 'ha': []},
                 'theilsen': {'ca': [], 'cb': [], 'c': [], 'n': [], 'ha': []},
                 'bayes': {'ca': [], 'cb': [], 'c': [], 'n': [], 'ha': []}}
    ca=[]
    ca_method=[]
    cb=[]
    cb_method=[]
    c=[]
    c_method=[]
    n=[]
    n_method=[]
    ha=[]
    ha_method=[]
    wrong_ref =[]
    tot=[]
    for jf in iter_json_files(args.directory, recursive=args.recursive, include_root=args.include_root):
        total += 1
        path, data, err = parse_json_file(jf)
        if err:
            errors += 1
            print(f"ERR | {path} | {err}")
            continue
        ok += 1
        #process_json(data, path)
        if len(data):
            entry_id = path.stem.split("_")[0]
            method = path.stem.split("_")[-1]
            if entry_id not in tot: tot.append(entry_id)
            for list_id in data:
                for atom in data[list_id]['offsets']:
                    if atom == 'ca' and abs(data[list_id]['offsets'][atom])>1.7:
                        if atom == 'ca' and entry_id not in wrong_ref: wrong_ref.append(entry_id)
                        print(f'{entry_id},{list_id},{method},{atom},{data[list_id]['offsets'][atom]}')
                    if atom == 'ca':
                        ca.append(data[list_id]['offsets'][atom])
                        ca_method.append(method)
                    elif atom == 'cb':
                        cb.append(data[list_id]['offsets'][atom])
                        cb_method.append(method)
                    elif atom == 'c':
                        c.append(data[list_id]['offsets'][atom])
                        c_method.append(method)
                    elif atom == 'n':
                        n.append(data[list_id]['offsets'][atom])
                        n_method.append(method)
                    elif atom == 'ha':
                        ha.append(data[list_id]['offsets'][atom])
                        ha_method.append(method)
                    dict_data[method][atom].append(data[list_id]['offsets'][atom])
    # df_ca = pd.DataFrame({'tukey':dict_data['tukey']['ca'],
    #                       'ransac':dict_data['ransac']['ca'],
    #                       'quantile':dict_data['quantile']['ca'],
    #                       'theilsen':dict_data['theilsen']['ca'],
    #                       'bayes':dict_data['bayes']['ca']})
    # df_cb = pd.DataFrame({'tukey':dict_data['tukey']['cb'],
    #                       'ransac':dict_data['ransac']['cb'],
    #                       'quantile':dict_data['quantile']['cb'],
    #                       'theilsen':dict_data['theilsen']['cb'],
    #                       'bayes':dict_data['bayes']['cb']})
    # df_c = pd.DataFrame({'tukey':dict_data['tukey']['c'],
    #                       'ransac':dict_data['ransac']['c'],
    #                       'quantile':dict_data['quantile']['c'],
    #                       'theilsen':dict_data['theilsen']['c'],
    #                       'bayes':dict_data['bayes']['c']})
    # df_n = pd.DataFrame({'tukey':dict_data['tukey']['n'],
    #                       'ransac':dict_data['ransac']['n'],
    #                       'quantile':dict_data['quantile']['n'],
    #                       'theilsen':dict_data['theilsen']['n'],
    #                       'bayes':dict_data['bayes']['n']})
    # df_ha = pd.DataFrame({'tukey':dict_data['tukey']['ha'],
    #                       'ransac':dict_data['ransac']['ha'],
    #                       'quantile':dict_data['quantile']['ha'],
    #                       'theilsen':dict_data['theilsen']['ha'],
    #                       'bayes':dict_data['bayes']['ha']})
    # fig = px.histogram(x=ca,color=ca_method,nbins=300,barmode='overlay',title='CA')
    # fig.show()
    # fig = px.histogram(x=cb,color=cb_method,barmode='overlay',title='CB')
    # fig.show()
    # fig = px.histogram(x=c,color=c_method,barmode='overlay',title='C')
    # fig.show()
    # fig = px.histogram(x=n,color=n_method,barmode='overlay',title='N')
    # fig.show()
    # fig = px.histogram(x=ha,color=ha_method,barmode='overlay',title='HA')
    # fig.show()

    print("-" * 80)
    print(f"Scanned: {total} file(s) | Parsed OK: {ok} | Errors: {errors}")
    for method in dict_data:
        for atom in dict_data[method]:
            print(f'{method},{atom},{len(dict_data[method][atom])},{np.mean(dict_data[method][atom])},{np.std(dict_data[method][atom])}')
    print (f'Wrong reference: {len(wrong_ref)}')
    print (f'Total: {len(tot)}')
    print (wrong_ref)

if __name__ == "__main__":
    main()
