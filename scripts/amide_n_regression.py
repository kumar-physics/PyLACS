#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Amide Nitrogen Regression Analysis
----------------------------------
Fits a robust linear model for amide nitrogen (N) chemical shifts using:
    y = ΔδN_i
    x = ΔδN_i − ΔδN_{i+1}
Pairs are formed only for consecutive residues (i, i+1) within the same
(Entity_assembly_ID, Entity_ID) "chain".

INPUT MODES
-----------
1) NMR-STAR (.str/.star) via pynmrstar (if installed)
   - Use --input-format star
   - ΔδN uses random-coil values from either:
       a) your uploaded random_coil.py (default, --rc-from module), or
       b) a CSV (--rc-from csv) via --rc-table with headers: resname,RC_N

2) CSV mode
   - Use --input-format csv
   - CSV must contain at least:
       entity_assembly_id, entity_id, comp_index_id, comp_id, N
     Optionally include RC_N; if present and selected, ΔδN = N - RC_N.

ROBUST METHODS
--------------
--method {ols,huber,theilsen,ransac}
Default is huber.

OUTPUTS
-------
- Prints slope, intercept (offset), and diagnostics.
- Optional: writes pairs to CSV (--export-csv) and a Plotly HTML scatter (--plot).

Author: ChatGPT (GPT-5 Thinking)
Date: 2025-09-12
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple, Any

# --- Optional imports (gracefully handled) ---
try:
    import pynmrstar  # type: ignore
except Exception:  # pragma: no cover
    pynmrstar = None

try:
    import numpy as np
except Exception as e:
    raise SystemExit("This script requires NumPy. Please install it: pip install numpy") from e

# scikit-learn for robust regressions
try:
    from sklearn.linear_model import LinearRegression, HuberRegressor, TheilSenRegressor, RANSACRegressor
except Exception as e:
    raise SystemExit("This script requires scikit-learn. Please install it: pip install scikit-learn") from e

# Plotly optional
try:
    import plotly.express as px
    import plotly.io as pio
except Exception:
    px = None
    pio = None

# Prefer the provided RandomCoil implementation
try:
    from random_coil import RandomCoil  # provided by user
except Exception:
    RandomCoil = None  # will error later if required


@dataclass(frozen=True)
class ResidueKey:
    entity_assembly_id: str
    entity_id: str
    comp_index_id: int
    comp_id: str  # three-letter code

    @property
    def chain_key(self) -> Tuple[str, str]:
        return (self.entity_assembly_id, self.entity_id)


@dataclass
class NitrogenDatum:
    key: ResidueKey
    N: float                 # observed amide nitrogen chemical shift
    rc_N: Optional[float]    # random-coil reference (if known)

    @property
    def delta_N(self) -> float:
        """ΔδN = observed - RC (if RC present), else treat observed as Δδ already."""
        if self.rc_N is not None:
            return float(self.N) - float(self.rc_N)
        return float(self.N)


def load_rc_table(rc_csv: Path) -> Dict[str, float]:
    """
    Load random-coil table CSV: expects columns 'resname','RC_N'
    Returns dict mapping 3-letter residue name (upper) -> RC_N float
    """
    rc: Dict[str, float] = {}
    with rc_csv.open("r", newline="") as f:
        for row in csv.DictReader(f):
            name = str(row.get("resname", "")).strip().upper()
            val = row.get("RC_N", None)
            if not name or val is None or str(val).strip() == "":
                continue
            try:
                rc[name] = float(val)
            except ValueError:
                continue
    if not rc:
        raise ValueError(f"No usable rows found in RC table: {rc_csv}")
    return rc


def normalize_tag(t: str) -> str:
    t = str(t).strip()
    if t.startswith("_"):
        t = t[1:]
    if "." in t:
        t = t.split(".")[-1]
    return t.lower()


def get_loop_tags(loop) -> List[str]:
    # Try several APIs to obtain tag names
    for attr in ("get_tag_names", "tags", "tag_names"):
        if hasattr(loop, attr):
            obj = getattr(loop, attr)
            names = obj() if callable(obj) else obj
            return [normalize_tag(n) for n in names]
    if hasattr(loop, "tag_order"):
        return [normalize_tag(n) for n in getattr(loop, "tag_order")]
    raise RuntimeError("Unable to obtain tag names from STAR loop")


def load_from_star(star_path: Path, rc_table: Optional[Dict[str, float]] = None) -> List[NitrogenDatum]:
    if pynmrstar is None:
        raise RuntimeError("pynmrstar is not available. Install it or use --input-format csv.")

    entry = pynmrstar.Entry.from_file(str(star_path))
    loops = entry.get_loops_by_category("Atom_chem_shift")

    records: List[NitrogenDatum] = []
    for loop in loops:
        cols = get_loop_tags(loop)
        name_to_idx = {name: i for i, name in enumerate(cols)}

        def first_present(cands: Sequence[str]) -> Optional[str]:
            for cand in cands:
                key = normalize_tag(cand)
                if key in name_to_idx:
                    return key
            return None

        col_atom_id   = first_present(["Atom_ID", "Atom_id", "Atom_ID_with_seq", "Atom_name"]) 
        col_val       = first_present(["Val", "Chem_shift_val", "Val_err", "Val_uncertainty", "Chemical_shift_value"]) 
        col_assem_id  = first_present(["Entity_assembly_ID", "Entity_assembly_id", "Assembly_ID", "Ensemble_ID"]) 
        col_ent_id    = first_present(["Entity_ID", "Entity_id"]) 
        col_comp_idx  = first_present(["Comp_index_ID", "Residue_seq_code", "Seq_ID", "Comp_index_id"]) 
        col_comp_id   = first_present(["Comp_ID", "Comp_id", "Residue_label", "Res_ID"]) 

        if not (col_atom_id and col_val and col_assem_id and col_ent_id and col_comp_idx and col_comp_id):
            continue  # skip non-conforming loops

        for row in loop.data:
            try:
                atom_name = str(row[name_to_idx[col_atom_id]]).strip().upper()
            except Exception:
                continue
            if atom_name not in {"N", "N15", "N15H"}:  # be tolerant to variants
                continue

            try:
                nval = float(str(row[name_to_idx[col_val]]))
            except Exception:
                continue

            try:
                assemb = str(row[name_to_idx[col_assem_id]]).strip()
                ent    = str(row[name_to_idx[col_ent_id]]).strip()
                comp_i = int(float(str(row[name_to_idx[col_comp_idx]]).strip()))
                comp   = str(row[name_to_idx[col_comp_id]]).strip().upper()
            except Exception:
                continue

            rcN = None
            if rc_table is not None:
                rcN = rc_table.get(comp.upper())

            key = ResidueKey(assemb, ent, comp_i, comp)
            records.append(NitrogenDatum(key=key, N=nval, rc_N=rcN))

    if not records:
        raise ValueError("No amide nitrogen (N) shifts found in the STAR file (category Atom_chem_shift).")
    return records


def load_from_csv(csv_path: Path,
                  col_map: Dict[str, str],
                  rc_in_csv: bool = True) -> List[NitrogenDatum]:
    """
    CSV must contain columns per col_map:
      entity_assembly_id, entity_id, comp_index_id, comp_id, N
    Optionally RC_N if rc_in_csv=True.
    """
    req = ["entity_assembly_id", "entity_id", "comp_index_id", "comp_id", "N"]
    for k in req:
        if k not in col_map:
            raise ValueError(f"Missing column mapping for '{k}'.")

    out: List[NitrogenDatum] = []
    with csv_path.open("r", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                assemb = str(row[col_map["entity_assembly_id"]]).strip()
                ent    = str(row[col_map["entity_id"]]).strip()
                comp_i = int(str(row[col_map["comp_index_id"]]).strip())
                comp   = str(row[col_map["comp_id"]]).strip().upper()
                nval   = float(str(row[col_map["N"]]).strip())
            except Exception:
                continue

            rcN = None
            if rc_in_csv and "RC_N" in col_map:
                rc_str = row.get(col_map["RC_N"], "")
                if rc_str not in (None, ""):
                    try:
                        rcN = float(str(rc_str))
                    except ValueError:
                        rcN = None

            key = ResidueKey(assemb, ent, comp_i, comp)
            out.append(NitrogenDatum(key=key, N=nval, rc_N=rcN))

    if not out:
        raise ValueError("No usable rows found in CSV.")
    return out


def build_neighbor_pairs(data: List[NitrogenDatum], rc_required: bool = False) -> Tuple[np.ndarray, np.ndarray, List[ResidueKey]]:
    """
    Build (x, y) pairs where:
        x = ΔδN_i − ΔδN_{i+1}
        y = ΔδN_i
    only when comp_index_id(i+1) exists in the same chain.
    """
    # Organize by chain
    from collections import defaultdict
    chains: Dict[Tuple[str, str], List[NitrogenDatum]] = defaultdict(list)
    for d in data:
        chains[d.key.chain_key].append(d)

    X: List[float] = []
    Y: List[float] = []
    tags: List[ResidueKey] = []

    for chain_key, residues in chains.items():
        residues_sorted = sorted(residues, key=lambda r: r.key.comp_index_id)
        for k in range(len(residues_sorted) - 1):
            i = residues_sorted[k]
            j = residues_sorted[k + 1]
            if j.key.comp_index_id != i.key.comp_index_id + 1:
                continue  # not consecutive

            if rc_required and (i.rc_N is None or j.rc_N is None):
                continue

            di = i.delta_N
            dj = j.delta_N
            X.append(di - dj)
            Y.append(di)
            tags.append(i.key)

    if len(X) == 0:
        raise ValueError("No consecutive (i, i+1) pairs found for nitrogens.")
    return np.asarray(X, dtype=float)[:, None], np.asarray(Y, dtype=float), tags


@dataclass
class FitResult:
    method: str
    slope: float
    intercept: float  # "offset"
    r2: Optional[float]
    inlier_mask: Optional[np.ndarray]
    n_pairs: int

    def as_dict(self) -> Dict[str, Any]:
        return {
            "method": self.method,
            "slope": self.slope,
            "intercept_offset": self.intercept,
            "r2": self.r2,
            "n_pairs": self.n_pairs,
            "n_inliers": int(self.inlier_mask.sum()) if self.inlier_mask is not None else self.n_pairs,
        }


def fit_model(X: np.ndarray, y: np.ndarray, method: str = "huber") -> FitResult:
    method = method.lower()
    n = len(y)

    if method == "ols":
        reg = LinearRegression()
        reg.fit(X, y)
        yhat = reg.predict(X)
        r2 = float(1.0 - np.sum((y - yhat) ** 2) / np.sum((y - y.mean()) ** 2)) if n > 1 else None
        return FitResult("ols", float(reg.coef_.ravel()[0]), float(reg.intercept_), r2, None, n)

    elif method == "huber":
        reg = HuberRegressor()
        reg.fit(X, y)
        yhat = reg.predict(X)
        # pseudo-R2
        ss_res = float(np.sum((y - yhat) ** 2))
        ss_tot = float(np.sum((y - y.mean()) ** 2)) if n > 1 else float("nan")
        r2 = float(1.0 - ss_res / ss_tot) if math.isfinite(ss_tot) and ss_tot > 0 else None
        return FitResult("huber", float(reg.coef_.ravel()[0]), float(reg.intercept_), r2, None, n)

    elif method == "theilsen":
        reg = TheilSenRegressor(random_state=0)
        reg.fit(X, y)
        yhat = reg.predict(X)
        ss_res = float(np.sum((y - yhat) ** 2))
        ss_tot = float(np.sum((y - y.mean()) ** 2)) if n > 1 else float("nan")
        r2 = float(1.0 - ss_res / ss_tot) if math.isfinite(ss_tot) and ss_tot > 0 else None
        return FitResult("theilsen", float(reg.coef_.ravel()[0]), float(reg.intercept_), r2, None, n)

    elif method == "ransac":
        base = LinearRegression()
        reg = RANSACRegressor(base_estimator=base, random_state=0)
        reg.fit(X, y)
        inlier_mask = reg.inlier_mask_
        est = reg.estimator_
        slope = float(est.coef_.ravel()[0])
        intercept = float(est.intercept_)
        # r2 on inliers only
        if inlier_mask is not None and inlier_mask.any():
            yhat = reg.predict(X[inlier_mask])
            y_in = y[inlier_mask]
            ss_res = float(np.sum((y_in - yhat) ** 2))
            ss_tot = float(np.sum((y_in - y_in.mean()) ** 2)) if y_in.size > 1 else float("nan")
            r2 = float(1.0 - ss_res / ss_tot) if math.isfinite(ss_tot) and ss_tot > 0 else None
        else:
            r2 = None
        return FitResult("ransac", slope, intercept, r2, inlier_mask, n)

    else:
        raise ValueError(f"Unknown method '{method}'. Choose from: ols, huber, theilsen, ransac")

def split_side_masks(X: np.ndarray):
    """Return boolean masks for positive and negative sides based on x = X[:,0]."""
    x = X[:, 0]
    m_pos = x >= 0.0
    m_neg = x < 0.0
    return m_pos, m_neg

def fit_both_sides(X: np.ndarray, y: np.ndarray, method: str):
    """Fit separate models for x>=0 and x<0; also fit a global model for reference.
    Returns (global_fit, pos_fit, neg_fit, pos_mask, neg_mask, pos_inlier_full, neg_inlier_full).
    """
    n = y.shape[0]
    m_pos, m_neg = split_side_masks(X)
    # Global reference fit
    global_fit = fit_model(X, y, method) if n >= 2 else None

    pos_fit = None
    neg_fit = None
    pos_inlier_full = None
    neg_inlier_full = None

    if m_pos.sum() >= 2:
        pos_fit = fit_model(X[m_pos], y[m_pos], method)
        pos_inlier_full = m_pos.copy()
        if pos_fit.inlier_mask is not None:
            mapped = pos_inlier_full.copy()
            mapped[:] = False
            mapped[m_pos] = pos_fit.inlier_mask
            pos_inlier_full = mapped
    if m_neg.sum() >= 2:
        neg_fit = fit_model(X[m_neg], y[m_neg], method)
        neg_inlier_full = m_neg.copy()
        if neg_fit.inlier_mask is not None:
            mapped = neg_inlier_full.copy()
            mapped[:] = False
            mapped[m_neg] = neg_fit.inlier_mask
            neg_inlier_full = mapped

    return global_fit, pos_fit, neg_fit, m_pos, m_neg, pos_inlier_full, neg_inlier_full





def save_pairs_csv(path: Path, X: np.ndarray, y: np.ndarray, tags: List[ResidueKey],
                   inlier_mask: Optional[np.ndarray] = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["entity_assembly_id", "entity_id", "comp_index_id", "comp_id", "x(ΔδN_i-ΔδN_{i+1})", "y(ΔδN_i)", "inlier"])
        for k, tag in enumerate(tags):
            inl = True if inlier_mask is None else bool(inlier_mask[k])
            w.writerow([tag.entity_assembly_id, tag.entity_id, tag.comp_index_id, tag.comp_id,
                        float(X[k, 0]), float(y[k]), int(inl)])


def make_plot_html(
    path: Path,
    X: np.ndarray,
    y: np.ndarray,
    fr_global: FitResult,
    fr_pos: Optional[FitResult],
    fr_neg: Optional[FitResult]
) -> None:
    if px is None:
        return
    import pandas as pd
    df = pd.DataFrame({"x": X[:, 0], "y": y})
    df["side"] = np.where(df["x"] >= 0, "pos", "neg")

    title = "Carbonyl regression (split sides)"
    fig = px.scatter(
        df, x="x", y="y", color="side", title=title,
        labels={"x": "ΔδCᵢ − ΔδCᵢ₊₁", "y": "ΔδCᵢ"}
    )

    if fr_pos is not None and (df["side"] == "pos").any():
        xseg = df.loc[df["side"] == "pos", "x"]
        xv = np.linspace(xseg.min(), xseg.max(), 100)
        yv = fr_pos.slope * xv + fr_pos.intercept
        fig.add_scatter(x=xv, y=yv, mode="lines",
                        name=f"pos fit: y={fr_pos.slope:.3f}x+{fr_pos.intercept:.3f}")

    if fr_neg is not None and (df["side"] == "neg").any():
        xseg = df.loc[df["side"] == "neg", "x"]
        xv = np.linspace(xseg.min(), xseg.max(), 100)
        yv = fr_neg.slope * xv + fr_neg.intercept
        fig.add_scatter(x=xv, y=yv, mode="lines",
                        name=f"neg fit: y={fr_neg.slope:.3f}x+{fr_neg.intercept:.3f}")

    path.parent.mkdir(parents=True, exist_ok=True)
    fig.write_html(str(path))


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Amide N robust regression: y=ΔδN_i vs x=ΔδN_i−ΔδN_{i+1}")
    p.add_argument("input", help="Input file: STAR (.str/.star) or CSV")
    p.add_argument("--input-format", choices=["star", "csv"], default=None,
                   help="Override autodetection of input format.")
    p.add_argument("--rc-table", type=str, default=None,
                   help="CSV with columns resname,RC_N (used for ΔδN = N - RC_N).")
    p.add_argument("--method", choices=["ols", "huber", "theilsen", "ransac"], default="huber")
    p.add_argument("--rc-policy", choices=["optional", "require", "none"], default="optional",
                   help="How to handle random-coil: optional (default), require (error/drop if missing), or none (treat observed N as ΔδN).")
    p.add_argument("--on-missing-rc", choices=["error", "drop"], default="error",
                   help="When --rc-policy=require: error if any RC is missing, or drop pairs that lack RC in i or i+1.")
    p.add_argument("--rc-from", choices=["module", "csv"], default="module",
                   help="Where to fetch random-coil values from. 'module': use uploaded random_coil.py; 'csv': use RC_N column / --rc-table.")
    p.add_argument("--rc-models", type=str, default="all",
                   help="Comma-separated model names for random_coil.RandomCoil (e.g., 'wis,wan,sch'). Use 'all' to average all models.")
    p.add_argument("--min-pairs", type=int, default=10, help="Minimum number of (i,i+1) pairs required.")
    p.add_argument("--outdir", type=str, default="amide_n_out")
    p.add_argument("--out", dest="outdir", type=str, help="Alias for --outdir")
    p.add_argument("--export-csv", action="store_true", help="Export x/y pairs and inlier mask to CSV.")
    p.add_argument("--plot", action="store_true", help="Write Plotly HTML scatter with fitted line.")
    # CSV column mapping
    p.add_argument("--csv-col-entity-assembly-id", default="entity_assembly_id")
    p.add_argument("--csv-col-entity-id", default="entity_id")
    p.add_argument("--csv-col-comp-index-id", default="comp_index_id")
    p.add_argument("--csv-col-comp-id", default="comp_id")
    p.add_argument("--csv-col-N", default="N")
    p.add_argument("--csv-col-RC_N", default="RC_N",
                   help="If present in CSV and you want to use it, keep name or map accordingly.")
    p.add_argument("--no-rc-in-csv", action="store_true",
                   help="If set, ignore RC_N column in CSV even if present.")
    return p.parse_args(argv)


def autodetect_format(path: Path) -> str:
    ext = path.suffix.lower()
    if ext in {".str", ".star"}:
        return "star"
    elif ext in {".csv", ".tsv"}:
        return "csv"
    return "csv"  # default to CSV if unknown


def main(argv: Optional[Sequence[str]] = None) -> None:
    args = parse_args(argv)
    in_path = Path(args.input).expanduser().resolve()
    if not in_path.exists():
        raise SystemExit(f"Input not found: {in_path}")

    fmt = args.input_format or autodetect_format(in_path)

    rc_table = None
    if args.rc_table:
        rc_table = load_rc_table(Path(args.rc_table))

    # Load data
    if fmt == "star":
        data = load_from_star(in_path, rc_table=rc_table if args.rc_from == "csv" else None)
    elif fmt == "csv":
        col_map = {
            "entity_assembly_id": args.csv_col_entity_assembly_id,
            "entity_id": args.csv_col_entity_id,
            "comp_index_id": args.csv_col_comp_index_id,
            "comp_id": args.csv_col_comp_id,
            "N": args.csv_col_N,
        }
        if not args.no_rc_in_csv and args.csv_col_RC_N:
            col_map["RC_N"] = args.csv_col_RC_N
            rc_in_csv = True
        else:
            rc_in_csv = False
        data = load_from_csv(in_path, col_map=col_map, rc_in_csv=rc_in_csv)
    else:
        raise SystemExit(f"Unsupported input format: {fmt}")

    # Build RandomCoil accessor if using module
    if args.rc_from == "module" and args.rc_policy != "none":
        if RandomCoil is None:
            raise SystemExit("random_coil.py not available, but --rc-from=module was requested.")
        rc_obj = RandomCoil()
        if args.rc_models.strip().lower() == "all":
            rc_models = None  # average across all models
        else:
            rc_models = [m.strip() for m in args.rc_models.split(",") if m.strip()]
        # attach RC to each datum using module
        for d in data:
            try:
                d.rc_N = rc_obj.get_value(d.key.comp_id, 'N', rc_models)
            except Exception:
                pass  # leave None; policy will handle it
    else:
        # If RC table provided (csv source), enrich missing RC_N values
        if rc_table is not None:
            for d in data:
                if d.rc_N is None:
                    d.rc_N = rc_table.get(d.key.comp_id.upper())

    # Enforce RC policy
    rc_required = (args.rc_policy == "require")
    if args.rc_policy == "none":
        for d in data:
            d.rc_N = None  # ensure delta_N uses observed N
    else:
        if rc_required:
            missing = sum(1 for d in data if d.rc_N is None)
            if missing > 0 and args.on_missing_rc == "error":
                raise SystemExit(f"Random-coil values missing for {missing} residues while --rc-policy=require. "
                                 f"Provide --rc-table or use --on-missing-rc=drop to drop such pairs.")

    # Build pairs
    X, y, tags = build_neighbor_pairs(data, rc_required=rc_required)
    if X.shape[0] < args.min_pairs:
        raise SystemExit(f"Only {X.shape[0]} pairs available (< min-pairs {args.min_pairs}).")

    # Split-side fits (and global reference)
    fr_global, fr_pos, fr_neg, m_pos, m_neg, inl_pos_full, inl_neg_full = fit_both_sides(X, y, args.method)

    # Report
    print("\nAmide Nitrogen Robust Regression (split sides)")
    print("--------------------------------")
    print(f"Total pairs   : {X.shape[0]} (pos={int(m_pos.sum())}, neg={int(m_neg.sum())})")
    if fr_global is not None:
        print(f"[global] slope={fr_global.slope:.6f}, intercept={fr_global.intercept:.6f}, R^2={fr_global.r2 if fr_global.r2 is not None else float('nan'):.6f}")
    if fr_pos is not None:
        print(f"[x>=0 ] slope={fr_pos.slope:.6f}, intercept={fr_pos.intercept:.6f}, R^2={fr_pos.r2 if fr_pos.r2 is not None else float('nan'):.6f}")
    else:
        print("[x>=0 ] insufficient points")
    if fr_neg is not None:
        print(f"[x<0  ] slope={fr_neg.slope:.6f}, intercept={fr_neg.intercept:.6f}, R^2={fr_neg.r2 if fr_neg.r2 is not None else float('nan'):.6f}")
    else:
        print("[x<0  ] insufficient points")

    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    # Save JSON summary
    summary = {
        "method": args.method,
        "counts": {"total_pairs": int(X.shape[0]), "pos_pairs": int(m_pos.sum()), "neg_pairs": int(m_neg.sum())},
        "fits": {
            "global": (fr_global.as_dict() if fr_global is not None else None),
            "pos": (fr_pos.as_dict() if fr_pos is not None else None),
            "neg": (fr_neg.as_dict() if fr_neg is not None else None)
        }
    }
    with (outdir / "summary.json").open("w") as f:
        json.dump(summary, f, indent=2)

    # Export pairs
    if args.export_csv:
        save_pairs_csv(outdir / "pairs.csv", X, y, tags, None, inl_pos_full, inl_neg_full)

    # Plot
    # Plot
    if args.plot:
        if px is None:
            print("Plotly not installed; skipping plot.", file=sys.stderr)
        else:
            make_plot_html(outdir / "scatter.html", X, y, fr_global, fr_pos, fr_neg)

    print(f"\nWrote outputs to: {outdir}")
    if args.export_csv:
        print(f"  - pairs.csv")
    print(f"  - summary.json")
    if args.plot and px is not None:
        print(f"  - scatter.html")


if __name__ == "__main__":
    main()
