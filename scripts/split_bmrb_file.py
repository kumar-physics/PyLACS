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




def split_entry(entry, input_name: str, outdir: Path, prefix: str = None) -> int:
    ent_tempate = copy.deepcopy(entry)
    for sf in entry:
        if sf.category== 'assigned_chemical_shifts':
            #print (sf.get_tag('ID'))
            ent_tempate.remove_saveframe(sf)
    count=0
    for sf in entry:
        if sf.category == 'assigned_chemical_shifts':
            new_ent = copy.deepcopy(ent_tempate)
            new_ent.add_saveframe(sf)
            list_id = sf.get_tag('ID')[0]
            ent_id = entry.entry_id
            stem = f'{ent_id}_{list_id}'
            if prefix:
                stem = f"{prefix}_{stem}"
            out_path = outdir / f"{stem}.str"

            # Write
            with open(out_path, "w", encoding="utf-8") as fh:
                fh.write(str(new_ent))

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
