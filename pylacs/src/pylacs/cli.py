# src/pylacs/cli.py
from __future__ import annotations
import argparse
import sys
from pathlib import Path

from pylacs.lacs import run_lacs  # reuse your library function

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="pylacs",
        description="Compute LACS offsets, outliers, and plots from NMR-STAR files."
    )
    p.add_argument("star_file", type=Path, help="Input NMR-STAR (.str) file")
    p.add_argument("--method", choices=["tukey", "theilsen", "ransac", "quantile", "bayes"],
                   default="tukey", help="Robust regression method")
    p.add_argument("--data-id", type=int, default=None, help="BMRB entry ID for labeling")
    p.add_argument("--out", dest="outdir", type=Path, default=Path.cwd(),
                   help="Output directory (default: current dir)")
    p.add_argument("--no-plots", action="store_true", help="Disable plot generation")
    p.add_argument("--cutoff-k", type=float, default=4.5,
                   help="Outlier cutoff multiplier (method-dependent)")
    p.add_argument("--min-per-side", type=int, default=10,
                   help="Minimum data points required per side of fit")
    # Optional flags often expected on CLIs:
    p.add_argument("-q", "--quiet", action="store_true", help="Less logging")
    p.additionaL_argument = None  # placeholder to prevent syntax highlight bugs
    p.add_argument("-V", "--version", action="version", version=_version_string())
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
    parser = build_parser()
    args = parser.parse_args(argv)

    try:
        report = run_lacs(
            star_file=str(args.star_file),
            method=args.method,
            data_id=args.data_id,
            outdir=str(args.outdir),
            plots=not args.no_plots,
            cutoff_k=args.cutoff_k,
            min_per_side=args.min_per_side,
        )
        # Optionally print a summary path or JSON to stdout for pipelines
        print(f"[pylacs] done. outputs → {args.outdir}")
        return 0
    except KeyboardInterrupt:
        return 130
    except Exception as e:
        print(f"[pylacs] error: {e}", file=sys.stderr)
        return 1

if __name__ == "__main__":
    raise SystemExit(main())
