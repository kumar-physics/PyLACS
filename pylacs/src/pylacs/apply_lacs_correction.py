#!/usr/bin/env python3
"""
apply_lacs_correction.py

Apply constant offsets to *selected* atoms among {CA, CB, C, N} in an
NMR-STAR file, restricted to a specific chemical-shift list ID, and append a
row to the Release/_Release loop documenting the correction.

Release behavior:
- Increments the release number (Ordinal / Release_number / Number) by +1.
- Copies Entry_ID from the previous Release row if present; else falls back to _Entry.ID.
- Writes a detailed offsets string (only for atoms actually applied) into the
  Release Detail/Details/Deatils column (prefers existing '_Release.Detail').
- Preserves existing loop schema; does not change tag order; pads unspecified cells with ".".
"""

import argparse
import sys
from typing import Dict, Optional, List, Tuple, Iterable
from datetime import date

import pynmrstar


# ----------------------------- Helpers -----------------------------

def _get_tag_names(loop) -> List[str]:
    if hasattr(loop, "tag_names") and loop.tag_names:
        return list(loop.tag_names)
    if hasattr(loop, "tags") and loop.tags:
        return list(loop.tags)
    if hasattr(loop, "columns") and loop.columns:
        return list(loop.columns)
    return []


def _get_loop_data(loop) -> List[List[str]]:
    if hasattr(loop, "data"):
        return loop.data
    if hasattr(loop, "rows"):
        return loop.rows
    return []


def _loop_is_pynmrstar_loop(obj) -> bool:
    return isinstance(obj, pynmrstar.Loop)


def _to_int_or_none(val) -> Optional[int]:
    try:
        if val is None:
            return None
        s = str(val).strip()
        if s in {"", ".", "?", "null", "None"}:
            return None
        return int(float(s))
    except Exception:
        return None


def _norm_atom_id(atom_id: str) -> str:
    return (atom_id or "").strip().upper()


def _norm_cat(cat: Optional[str]) -> str:
    return (cat or "").strip().lstrip("_").lower()


# ------------------------- List ID resolution -------------------------

def resolve_list_id_for_loop(sf: pynmrstar.Saveframe, loop) -> Optional[int]:
    if not _loop_is_pynmrstar_loop(loop):
        return None

    tag_names = _get_tag_names(loop)
    data = _get_loop_data(loop)

    loop_id_candidates = [
        "Assigned_chem_shift_list_ID",
        "Chem_shift_list_ID",
        "List_ID",
        "assigned_chem_shift_list_ID",
        "chem_shift_list_ID",
        "list_ID",
    ]
    for col in loop_id_candidates:
        if col in tag_names:
            idx = tag_names.index(col)
            if data:
                return _to_int_or_none(data[0][idx])

    sf_tag_candidates = [
        "Assigned_chem_shift_list.ID",
        "Chem_shift_list.ID",
        "Assigned_chem_shift_list_ID",
        "Chem_shift_list_ID",
        "ID",
    ]
    for tag in sf_tag_candidates:
        try:
            return _to_int_or_none(sf.get_tag(tag))
        except KeyError:
            continue

    return None


# ------------------------------ Core correction ------------------------------

def apply_offsets_to_entry_for_list(entry: pynmrstar.Entry,
                                    offsets: Dict[str, float],
                                    target_list_id: int,
                                    atoms_to_apply: Iterable[str]) -> Dict[str, int]:
    """
    Apply offsets to Atom_chem_shift loops in-place, but only for the specified
    chemical-shift list ID and only for atom names in `atoms_to_apply` (subset of CA, CB, C, N).
    """
    atoms_set = {a.upper() for a in atoms_to_apply}
    # Canonical set we support
    allowed = {"CA", "CB", "C", "N"}

    counts = {k: 0 for k in allowed}
    counts.update({"total": 0, "loops_seen": 0, "loops_matched": 0})

    def eligible_atom(atom: str) -> Optional[str]:
        a = _norm_atom_id(atom)
        if a in allowed and a in atoms_set:
            return a
        return None

    for sf in entry:
        cat = getattr(sf, "category", "") or ""
        if _norm_cat(cat) not in {"assigned_chemical_shifts", "chemical_shifts"}:
            continue

        try:
            loops = sf.get_loops_by_category("Atom_chem_shift")
        except Exception:
            loops = getattr(sf, "loops", [])

        for loop in loops or []:
            if not _loop_is_pynmrstar_loop(loop):
                continue

            counts["loops_seen"] += 1

            list_id = resolve_list_id_for_loop(sf, loop)
            if list_id is None or list_id != target_list_id:
                continue
            counts["loops_matched"] += 1

            tags = _get_tag_names(loop)
            rows = _get_loop_data(loop)

            atom_col = None
            for name in ("Atom_ID", "Atom_name", "Atom_name_remap", "Atom_id", "atom_ID", "atom_name"):
                if name in tags:
                    atom_col = tags.index(name)
                    break

            val_col = None
            for name in ("Val", "Chem_shift_val", "Chem_shift_value", "val", "Chem_shift_value_val"):
                if name in tags:
                    val_col = tags.index(name)
                    break

            if atom_col is None or val_col is None:
                continue

            for row in rows:
                key = eligible_atom(row[atom_col])
                if key is None:
                    continue
                off = float(offsets.get(key, 0.0))
                if off == 0.0:
                    continue
                try:
                    old_val = float(row[val_col])
                except (TypeError, ValueError):
                    continue
                row[val_col] = f"{(old_val + off):.6g}"
                counts[key] += 1
                counts["total"] += 1

    return counts


# ----------------------------- Release helpers -----------------------------

def _find_entry_information_sf(entry: pynmrstar.Entry) -> Optional[pynmrstar.Saveframe]:
    try:
        eifs = entry.get_saveframes_by_category("entry_information")
        if eifs:
            return eifs[0]
    except Exception:
        pass
    for sf in entry:
        if _norm_cat(getattr(sf, "category", "")) == "entry_information":
            return sf
    for sf in entry:
        return sf
    return None


def _find_release_loop(sf: pynmrstar.Saveframe) -> Optional[pynmrstar.Loop]:
    for cat in ("Release", "_Release"):
        try:
            loops = sf.get_loops_by_category(cat)
            if loops:
                return loops[0]
        except Exception:
            pass
    for lp in (getattr(sf, "loops", []) or []):
        if _loop_is_pynmrstar_loop(lp) and _norm_cat(getattr(lp, "category", "")) == "release":
            return lp
    return None


def _create_release_loop(sf: pynmrstar.Saveframe) -> pynmrstar.Loop:
    existing = _find_release_loop(sf)
    if existing is not None:
        return existing
    # Minimal, underscored style to match many BMRB files
    try:
        loop = pynmrstar.Loop.from_scratch("_Release", ["_Release.Author", "_Release.Date", "_Release.Detail"])
    except TypeError:
        loop = pynmrstar.Loop()
        try:
            loop.category = "_Release"
        except Exception:
            pass
        for t in ["_Release.Author", "_Release.Date", "_Release.Detail"]:
            if hasattr(loop, "add_tag"):
                loop.add_tag(t)
            else:
                if hasattr(loop, "tags"):
                    loop.tags.append(t)
    sf.add_loop(loop)
    return loop


def _first_present(names: List[str], tags: List[str]) -> Optional[int]:
    for n in names:
        if n in tags:
            return tags.index(n)
    return None


def _read_entry_id_from_entry_info(entry: pynmrstar.Entry) -> Optional[str]:
    sf = _find_entry_information_sf(entry)
    if sf is None:
        return None
    candidates = ["_Entry.ID", "Entry.ID", "ID", "_Entry.Sf_framecode", "Entry.Sf_framecode"]
    for t in candidates:
        try:
            v = sf.get_tag(t)
            s = str(v).strip()
            if s and s not in {".", "?"}:
                return s
        except KeyError:
            continue
    return None


# ----------------------------- Release append -----------------------------

def append_release_note(entry: pynmrstar.Entry,
                        author: str,
                        details: str,
                        date_str: Optional[str] = None) -> None:
    """
    Append a row to the existing Release/_Release loop (or create if missing).
    - Keeps existing schema/columns.
    - Increments release number.
    - Copies Entry_ID from previous row; else falls back to Entry.ID.
    - Writes `details` string into the Detail/Details/Deatils field (prefers 'Detail' if present).
    """
    date_str = date_str or date.today().isoformat()
    sf = _find_entry_information_sf(entry)
    if sf is None:
        return

    release_loop = _find_release_loop(sf)
    if release_loop is None:
        release_loop = _create_release_loop(sf)

    tags = _get_tag_names(release_loop)
    rows = _get_loop_data(release_loop)

    # Build a full-width row of '.'
    row = ["." for _ in tags]

    # Column indices (prefer underscored names if present)
    author_idx  = _first_present(["_Release.Author", "Release.Author", "Author"], tags)
    date_idx    = _first_present(["_Release.Date", "Release.Date", "Date"], tags)

    # Prefer singular 'Detail'; fall back to Details/Deatils
    detail_idx = _first_present(
        ["_Release.Detail", "Release.Detail", "Detail",
         "_Release.Details", "Release.Details", "Details",
         "_Release.Deatils", "Release.Deatils", "Deatils"],
        tags
    )

    # Fill core fields if present
    if author_idx is not None:
        row[author_idx] = author
    if date_idx is not None:
        row[date_idx] = date_str
    if detail_idx is not None:
        row[detail_idx] = details  # write the details string here

    # Increment release number if a numeric column exists
    release_num_idx = _first_present(
        ["_Release.Ordinal", "Release.Ordinal", "Ordinal",
         "_Release.Release_number", "Release.Release_number", "Release_number",
         "_Release.Number", "Release.Number", "Number"],
        tags
    )
    if release_num_idx is not None:
        max_num = 0
        for r in rows:
            try:
                v = str(r[release_num_idx]).strip()
                if v not in {"", ".", "?"}:
                    max_num = max(max_num, int(float(v)))
            except Exception:
                continue
        row[release_num_idx] = str(max_num + 1)

    # Entry_ID: copy from previous, else entry_information fallback
    entry_id_idx = _first_present(["_Release.Entry_ID", "Release.Entry_ID", "Entry_ID"], tags)
    if entry_id_idx is not None:
        prev_entry_id = None
        for r in reversed(rows):
            val = str(r[entry_id_idx]).strip()
            if val and val not in {".", "?"}:
                prev_entry_id = val
                break
        if prev_entry_id is None:
            prev_entry_id = _read_entry_id_from_entry_info(entry) or "."
        row[entry_id_idx] = prev_entry_id

    # Optional: set Type='Correction' if column exists
    type_idx = _first_present(["_Release.Type", "Release.Type", "Type"], tags)
    if type_idx is not None:
        row[type_idx] = "Correction"

    # Append row with exact width
    if hasattr(release_loop, "add_data"):
        release_loop.add_data([row])
    else:
        _get_loop_data(release_loop).append(row)


# --------------------------- Public wrapper API ---------------------------

def apply_selected_offsets_and_note(
    input_path: str,
    output_path: str,
    list_id: int,
    offsets: Dict[str, float],
    atoms: Iterable[str] = ("CA", "CB", "C", "N"),
    release_author: str = "BMRB",
    release_details: Optional[str] = None,
) -> Dict[str, int]:
    """
    Single-call wrapper for external use (e.g., from lacs.py or notebooks).

    Parameters
    ----------
    input_path : str
        Path to input NMR-STAR file.
    output_path : str
        Where to write corrected file.
    list_id : int
        Chemical-shift list ID to modify.
    offsets : Dict[str, float]
        Offsets per atom; only keys in {'CA','CB','C','N'} are used.
    atoms : Iterable[str]
        Subset of atoms to which offsets should be applied (defaults to all four).
    release_author : str
        Author to place in the Release row (default 'BMRB').
    release_details : Optional[str]
        If provided, used verbatim in the Release Detail field; otherwise composed automatically.

    Returns
    -------
    Dict[str, int]
        Counts dictionary with per-atom and total updates.
    """
    # Normalize offsets and atom list
    allowed = {"CA", "CB", "C", "N"}
    atoms_upper = [a.upper() for a in atoms if a is not None]
    atoms_use = [a for a in atoms_upper if a in allowed]
    offs = {k: float(offsets.get(k, 0.0)) for k in allowed}
    print (input_path)
    entry = pynmrstar.Entry.from_file(input_path)
    counts = apply_offsets_to_entry_for_list(entry, offs, list_id, atoms_use)

    # Compose default details message if not provided
    if release_details is None:
        parts = [f"{a}={offs[a]:+g}" for a in atoms_use]
        detail_offsets = ", ".join(parts) if parts else "none"
        details = (
            f"Applied chemical-shift offsets (ppm) to list_id {list_id}: "
            f"{detail_offsets}. Rows updated: " +
            ", ".join(f"{a}={counts[a]}" for a in atoms_use) +
            f"; total={counts['total']}."
        )
    else:
        details = release_details

    append_release_note(entry, author=release_author, details=details)
    entry.write_to_file(output_path)
    return counts


# ------------------------------------ CLI --------------------------------------

def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Apply selected CA/CB/C/N offsets to a specific chemical-shift list in an NMR-STAR file and record a Release note."
    )
    p.add_argument("input", help="Input NMR-STAR file (e.g., .str)")
    p.add_argument("output", help="Output NMR-STAR file (corrected)")
    p.add_argument("--list-id", required=True, type=int,
                   help="Chemical-shift list ID to modify (integer)")
    p.add_argument("--ca", type=float, default=0.0, help="Offset to add to CA shifts (ppm)")
    p.add_argument("--cb", type=float, default=0.0, help="Offset to add to CB shifts (ppm)")
    p.add_argument("--c",  type=float, default=0.0, help="Offset to add to C (carbonyl) shifts (ppm)")
    p.add_argument("--n",  type=float, default=0.0, help="Offset to add to N shifts (ppm)")

    # NEW: restrict which atoms to apply
    p.add_argument(
        "--atoms",
        nargs="+",
        choices=["CA", "CB", "C", "N", "ca", "cb", "c", "n"],
        default=["CA", "CB", "C", "N"],
        help="Subset of atoms to correct (choose from CA CB C N). Default: all four."
    )

    p.add_argument("--dry-run", action="store_true",
                   help="Report intended changes but do not write output.")
    return p.parse_args()


def main():
    args = _parse_args()
    # Normalize selected atoms to upper
    atoms = [a.upper() for a in args.atoms]

    offsets = {"CA": args.ca, "CB": args.cb, "C": args.c, "N": args.n}

    if all(abs(float(offsets[k])) < 1e-15 for k in offsets):
        print("All offsets are zero; no changes will be applied.", file=sys.stderr)

    try:
        entry = pynmrstar.Entry.from_file(args.input)
    except Exception as e:
        print(f"ERROR: Failed to read input NMR-STAR file '{args.input}': {e}", file=sys.stderr)
        sys.exit(1)

    counts = apply_offsets_to_entry_for_list(entry, offsets, args.list_id, atoms)

    # Compose the text that goes into Release.Detail (only selected atoms)
    parts = [f"{a}={offsets[a]:+g}" for a in atoms]
    detail_offsets = ", ".join(parts) if parts else "none"
    details_text = (
        f"Applied chemical-shift offsets (ppm) to list_id {args.list_id}: "
        f"{detail_offsets}. Rows updated: " +
        ", ".join(f"{a}={counts[a]}" for a in atoms) +
        f"; total={counts['total']}."
    )
    append_release_note(entry, author="BMRB", details=details_text)

    # Report
    print("Target chemical-shift list ID:", args.list_id)
    print("Atoms corrected:", ", ".join(atoms) if atoms else "(none)")
    print("Offset summary (ppm):")
    for a in atoms:
        print(f"  {a}: {float(offsets[a]):+g}")
    print("Rows updated:")
    for a in atoms:
        print(f"  {a}: {counts[a]}")
    print(f"Total updated: {counts['total']}")
    print(f"Loops scanned: {counts['loops_seen']}, matched list ID: {counts['loops_matched']}")

    if args.dry_run:
        print("\nDry run requested; no file written.")
        return

    try:
        entry.write_to_file(args.output)
        print(f"Corrected NMR-STAR written to: {args.output}")
    except Exception as e:
        print(f"ERROR: Failed to write output file '{args.output}': {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
