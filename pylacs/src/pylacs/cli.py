# src/pylacs/cli.py
from __future__ import annotations
import argparse
import sys
from pathlib import Path
import json

from pylacs.lacs import run_lacs  # reuse your library function
from pylacs.apply_lacs_correction import apply_selected_offsets_and_note
from pylacs.lacs import _extract_offsets_for_list,_default_json_path,_normalize_atoms

def build_parser() -> argparse.ArgumentParser:
    """
    Build a command-line argument parser for pylacs.

    :param argv: List of command-line arguments (default: sys.argv[1:])
    :return: ArgumentParser object
    """
    p = argparse.ArgumentParser(
        prog="pylacs",
        description="Compute LACS offsets, outliers, and plots from NMR-STAR files."
    )
    p.add_argument("star_file", help="Path to NMR-STAR .str file")
    p.add_argument("--method", default='bayes',
                   choices=["tukey", "theilsen", "ransac", "quantile", "bayes"],
                   help="Robust regression method to use (default Bayes)")
    p.add_argument("--data-id", default="LACS")
    p.add_argument("--rc-model", nargs="*", default=None,
                   help="Random-coil model alias(es), e.g. wis wan; omit for average of all")
    p.add_argument("--out", type=Path, default=None, help="Output directory for plots")
    p.add_argument("--no-plots", action="store_true", help="Disable plotting entirely")
    p.add_argument('--min-per-side', type=int, default=5,
                   help='Minimum number of points required on each sign side (default: 5).')
    p.add_argument("--cutoff-k", type=float, default=5.0, help="Outlier cutoff multiplier (default 5.0)")
    p.add_argument("--apply-offsets", action="store_true",
                   help="After computing offsets, apply them to the input .str.")
    p.add_argument("--out-format", choices=["json", "star", "both"], default="both",
                   help="Output format (default: both)")
    p.add_argument("--json-out", default=None, help="Explicit JSON path (overrides default).")
    p.add_argument("--star-out", default=None, help="Explicit STAR path (overrides default).")
    p.add_argument("--apply-corrections", action="store_true",
                   help="Apply LACS offsets to input STAR and write a corrected file.")
    p.add_argument("--correction-atoms", nargs="+", default=["CA", "CB", "C", "N"],
                   help="Atoms to correct when applying offsets.")
    p.add_argument("--release-author", default="BMRB",
                   help="Author to record in Release loop for corrections.")
    p.add_argument("--output-corrected", default=None,
                   help="Path to corrected output STAR (.str) file.")
    return p

def _version_string() -> str:
    try:
        from importlib.metadata import version
        v = version("pylacs")
    except Exception:
        # fallback to package attribute if editable install
        try:
            from pylacs import __version__ as v
        except Exception:
            v = "unknown"
    return f"pylacs {v}"

def main(argv: list[str] | None = None) -> int:
    """Main entry point for pylacs command-line interface."""
    parser = build_parser()
    args = parser.parse_args(argv)
    rc_model = None if args.rc_model == [] else (args.rc_model if args.rc_model is not None else None)

    try:
        report = run_lacs(args.star_file, method=args.method, data_id=args.data_id,
                          rc_model=rc_model, outdir=args.out, plots=not args.no_plots, cutoff_k=args.cutoff_k,
                          min_per_side=args.min_per_side,
                          write_format=args.out_format,
                          json_out=Path(args.json_out) if args.json_out else None,
                          star_out=Path(args.star_out) if args.star_out else None,
                          params_for_star=vars(args),  # capture CLI args as metadata
                          apply_corrections=bool(args.apply_corrections),
                          correction_atoms=args.correction_atoms,
                          release_author=args.release_author,
                          output_corrected=(Path(args.output_corrected) if args.output_corrected else None),
                          )
    #     # Optionally print a summary path or JSON to stdout for pipelines
    #     print(f"[pylacs] done. outputs → {args.outdir}")
    #     json_path = args.json_out or _default_json_path(args.outdir, args.data_id, args.method)
    #     json_path.parent.mkdir(parents=True, exist_ok=True)
    #     with open(json_path, "w", encoding="utf-8") as f:
    #         json.dump(report, f, indent=2)
    #     if args.apply_offsets:
    #         if apply_selected_offsets_and_note is None:
    #             raise SystemExit("Cannot apply offsets: apply_lacs_correction.py not importable.")
    #         if not args.output_corrected:
    #             corrected_path = json_path.parent / f'{args.data_id}_corrected,str'
    #             corrected_path.parent.mkdir(parents=True, exist_ok=True)
    #             #raise SystemExit("--output-corrected is required with --apply-offsets.")
    #         else:
    #             corrected_path = args.output_corrected
    #         # print (report)
    #         for list_id in report:
    #             try:
    #                 offsets_uc = _extract_offsets_for_list(report, list_id)
    #             except KeyError as e:
    #                 # Graceful fallback: try reading the freshly written JSON (identical content)
    #                 with open(json_path, "r", encoding="utf-8") as f:
    #                     report2 = json.load(f)
    #                 offsets_uc = _extract_offsets_for_list(report2, list_id)
    #
    #             atoms_use = _normalize_atoms(args.atoms)
    #             parts = [f"{a}={offsets_uc[a]:+g}" for a in atoms_use]
    #             details = (
    #                     f"LACS correction applied to list_id {list_id} ({args.data_id}): "
    #                     + (", ".join(parts) if parts else "none")
    #                     + ". Source: LACS (same run)."
    #             )
    #             if len(report) > 1:
    #                 list_ids = list(report.keys())
    #                 if list_ids.index(list_id) > 0:
    #                     counts = apply_selected_offsets_and_note(
    #                         input_path=str(args.star_file),
    #                         output_path=corrected_path,
    #                         list_id=int(list_id),
    #                         offsets=offsets_uc,
    #                         atoms=atoms_use,
    #                         release_author=args.release_author,
    #                         release_details=details,
    #                     )
    #                 else:
    #                     counts = apply_selected_offsets_and_note(
    #                         input_path=str(args.star_file),
    #                         output_path=corrected_path,
    #                         list_id=int(list_id),
    #                         offsets=offsets_uc,
    #                         atoms=atoms_use,
    #                         release_author=args.release_author,
    #                         release_details=details,
    #                     )
    #             else:
    #                 counts = apply_selected_offsets_and_note(
    #                     input_path=str(args.star_file),
    #                     output_path=corrected_path,
    #                     list_id=int(list_id),
    #                     offsets=offsets_uc,
    #                     atoms=atoms_use,
    #                     release_author=args.release_author,
    #                     release_details=details,
    #                 )
    #         print("Applied offsets to file:")
    #         for a in atoms_use:
    #             print(f"  {a}: {counts[a]}")
    #         print(f"Total updated: {counts['total']}")
    #         print(f"Wrote corrected file: {args.output_corrected}")
    #     return 0
    except KeyboardInterrupt:
        return 130
    except Exception as e:
        print(f"[pylacs] error: {e}", file=sys.stderr)
        return 1

if __name__ == "__main__":
    raise SystemExit(main())
