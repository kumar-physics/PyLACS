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
    print(f"OK  | {path} | keys: [{preview}]")


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

    for jf in iter_json_files(args.directory, recursive=args.recursive, include_root=args.include_root):
        total += 1
        path, data, err = parse_json_file(jf)
        if err:
            errors += 1
            print(f"ERR | {path} | {err}")
            continue
        ok += 1
        process_json(data, path)

    print("-" * 80)
    print(f"Scanned: {total} file(s) | Parsed OK: {ok} | Errors: {errors}")


if __name__ == "__main__":
    main()
