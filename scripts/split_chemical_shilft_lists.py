#!/usr/bin/env python3
"""
split_chemical_shift_lists.py

Split an NMR-STAR file that contains multiple chemical shift lists
(Assigned_chem_shift_list saveframes) into separate .str files.

Each output retains ALL other saveframes/metadata; only one Assigned_chem_shift_list
saveframe is kept per output file.

Usage:
  python split_chemical_shift_lists.py input.str -o outdir
"""

import argparse
import copy
import logging
from pathlib import Path

import pynmrstar


def _all_saveframes(entry):
    """Return all saveframes using the public API (no .saveframes attr)."""
    try:
        return list(entry.get_saveframes())
    except Exception:
        # Fallback for very old versions (unlikely, but safe)
        return list(getattr(entry, "frame_list", []))


def find_assigned_chem_shift_sfs(entry):
    """
    Find 'Assigned_chem_shift_list' saveframes (STAR v3).
    Be tolerant to naming differences / older STAR variants.
    """
    # Preferred: category-based lookup
    for cat in ("assigned_chem_shift_list", "assigned_chemical_shifts"):
        try:
            sfs = entry.get_saveframes_by_category(cat)
            if sfs:
                return list(sfs)
        except Exception:
            pass

    # Fallback: scan all saveframes and match category or name heuristics
    hits = []
    for sf in _all_saveframes(entry):
        cat = (getattr(sf, "category", "") or "").lower()
        name = (getattr(sf, "name", "") or "").lower()
        if cat in ("assigned_chem_shift_list", "assigned_chemical_shifts"):
            hits.append(sf)
            continue
        if name.startswith("assigned_chem_shift_list") or name.startswith("assigned_chemical_shifts"):
            hits.append(sf)
    return hits


def sanitize_for_filename(s: str) -> str:
    bad = '<>:"/\\|?* \t\r\n'
    return "".join('_' if ch in bad else ch for ch in s)


def derive_entry_id(entry, fallback: str) -> str:
    # Try standard entry_id
    entry_id = getattr(entry, "entry_id", None)
    if entry_id:
        return str(entry_id)
    # Try datablock name
    for attr in ("get_datablock_name", "get_data_block_name"):
        if hasattr(entry, attr):
            try:
                dbn = getattr(entry, attr)()
                if dbn:
                    return dbn.replace("data_", "")
            except Exception:
                pass
    # Fallback to input filename stem
    return Path(fallback).stem


def split_entry(entry, input_name: str, outdir: Path, prefix: str = None) -> int:
    shift_sfs = find_assigned_chem_shift_sfs(entry)
    if not shift_sfs:
        logging.warning("No Assigned_chem_shift_list saveframes found. Nothing to split.")
        return 0

    entry_id = derive_entry_id(entry, input_name)
    count = 0

    for idx, keep_sf in enumerate(shift_sfs, start=1):
        # Deep copy whole entry to retain all metadata
        new_entry = copy.deepcopy(entry)

        # Collect all shift-list frames in the copy
        to_remove = []
        for sf in _all_saveframes(new_entry):
            cat = (getattr(sf, "category", "") or "").lower()
            name = (getattr(sf, "name", "") or "")
            is_shift = (
                cat in ("assigned_chem_shift_list", "assigned_chemical_shifts")
                or name.lower().startswith("assigned_chem_shift_list")
                or name.lower().startswith("assigned_chemical_shifts")
            )
            if is_shift and name != getattr(keep_sf, "name", None):
                to_remove.append(sf)

        # Remove all other shift lists
        for sf in to_remove:
            try:
                new_entry.remove_saveframe(sf)
            except Exception:
                # Older versions also support deleting by object; if needed, try by name
                try:
                    new_entry.remove_saveframe(getattr(sf, "name", ""))
                except Exception:
                    pass  # keep going—worst case an extra shift list remains

        # Build output path
        keep_name = getattr(keep_sf, "name", f"shiftlist_{idx}") or f"shiftlist_{idx}"
        stem = f"{entry_id}__{sanitize_for_filename(keep_name)}"
        if prefix:
            stem = f"{prefix}__{stem}"
        out_path = outdir / f"{stem}.str"

        # Write
        with open(out_path, "w", encoding="utf-8") as fh:
            fh.write(str(new_entry))

        logging.info("Wrote %s", out_path)
        count += 1

    return count


def main():
    ap = argparse.ArgumentParser(description="Split multiple chemical shift lists into separate NMR-STAR files.")
    ap.add_argument("input", help="Input NMR-STAR .str file (BMRB entry).")
    ap.add_argument("-o", "--outdir", required=True, help="Output directory for split .str files.")
    ap.add_argument("--prefix", help="Optional prefix for output filenames.", default=None)
    ap.add_argument("-q", "--quiet", action="store_true", help="Reduce logging.")
    args = ap.parse_args()

    logging.basicConfig(
        level=logging.WARNING if args.quiet else logging.INFO,
        format="%(levelname)s: %(message)s"
    )

    in_path = Path(args.input)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    try:
        entry = pynmrstar.Entry.from_file(str(in_path))
    except Exception as e:
        logging.error("Failed to read input STAR file: %s", e)
        raise SystemExit(1)

    n = split_entry(entry, input_name=in_path.name, outdir=outdir, prefix=args.prefix)
    if n == 0:
        logging.warning("No output generated.")
    else:
        logging.info("Done. Generated %d file(s).", n)


if __name__ == "__main__":
    main()
