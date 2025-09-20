#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
LACS unified CLI (H/N variant)
==============================

This version computes robust linear fits for Δδ(H) and Δδ(N) using

    x = ΔδN − ΔδH

For each of the two nuclei (H and N), the data are split by the sign of x
(x ≥ 0 vs x < 0), a line is fit on each side, and the **offset** is taken as
the negative average of the two intercepts (or the single available intercept
if only one side has enough points). Outliers are flagged via robust residual
scoring with a smooth probability map.

Supported methods (unchanged):
- ``tukey``: Tukey biweight RLM (statsmodels)
- ``theilsen``: Theil–Sen median-of-slopes (scikit-learn)
- ``ransac``: RANSAC with MAD-based threshold (scikit-learn)
- ``quantile``: Quantile regression at τ=0.5 (statsmodels)
- ``bayes``: Bayesian Student-t regression (PyMC) with offset uncertainties

Usage
-----
.. code-block:: bash

    python lacs_unified_hn.py ENTRY.str --method tukey --data-id myprot --out figs --json-out myprot_tukey.json

Dependencies
------------
- *Common:* ``pynmrstar``, your ``random_coil`` provider on PYTHONPATH
- *Plotting (optional):* ``plotly`` (and ``kaleido`` for PDF export)
- *By method:*
  - Tukey/Quantile: ``statsmodels``
  - Theil–Sen/RANSAC: ``scikit-learn``
  - Bayesian: ``pymc`` (PyMC)
"""
from __future__ import annotations

import json
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

# Minimum points required per sign side
MIN_PER_SIDE_DEFAULT = 5

# ----------------------------- Optional plotting ------------------------------
try:  # pragma: no cover - optional
    import plotly.express as px
    import plotly.graph_objects as go
    PLOTLY_AVAILABLE = True
except Exception:  # pragma: no cover - optional
    PLOTLY_AVAILABLE = False

# ----------------------------- External deps ----------------------------------
try:
    import pynmrstar
except Exception as e:
    raise SystemExit("This script requires 'pynmrstar' (pip install pynmrstar).") from e

try:
    from pylacs.random_coil import RandomCoil
except Exception as e:
    raise SystemExit("Couldn't import random_coil.py. Ensure it is on PYTHONPATH or in the same folder.") from e

# The correction utility still supports only CA/CB/C/N as before.
try:
    from pylacs.apply_lacs_correction import apply_selected_offsets_and_note
except Exception as e:
    apply_selected_offsets_and_note = None  # allow running fits/plots/JSON without apply
try:
    # Also allow a flat import if present in cwd
    from apply_lacs_correction import apply_selected_offsets_and_note as _alt_apply
    apply_selected_offsets_and_note = apply_selected_offsets_and_note or _alt_apply
except Exception:
    pass


# ----------------------- helpers for offset application -----------------------

def _extract_offsets_for_list(report: Dict, list_id: int) -> Dict[str, float]:
    """
    Given a full LACS report (dict) extract offsets for the specified list_id.
    Returns keys in UPPER CASE: {'CA','CB','C','N'} (missing → 0.0).
    Note: H offsets are NOT applied by the current STAR correction tool.
    """
    sid = str(list_id)
    if sid not in report:
        # if the caller passed only the inner object (one list), accept that too
        candidate = report.get("offsets")
        if isinstance(candidate, dict):
            offs = candidate
        else:
            raise KeyError(f"List id {list_id} not found in report keys {list(report.keys())}")
    else:
        offs = report[sid].get("offsets", {})

    out = {"CA": 0.0, "CB": 0.0, "C": 0.0, "N": 0.0}
    for k_lc, v in offs.items():
        k_up = k_lc.upper()
        if k_up in out:
            try:
                out[k_up] = float(v)
            except Exception:
                pass
    return out


def _normalize_atoms(atoms: Sequence[str] | None) -> List[str]:
    # keep only the atoms supported by the STAR correction routine
    allowed = {"CA", "CB", "C", "N"}
    if not atoms:
        return ["CA", "CB", "C", "N"]
    atoms_up = [a.upper() for a in atoms]
    return [a for a in atoms_up if a in allowed]


# ---------------------------------------------------------------------------
# PUBLIC WRAPPER: apply correction from an existing JSON report (unchanged)
# ---------------------------------------------------------------------------

def apply_offset_correction(
    str_file: str,
    data_id: str,
    lacs_output: str,
    output_corrected: str,
    list_id: int,
    atoms: Sequence[str] = ("CA", "CB", "C", "N"),
    release_author: str = "BMRB",
) -> Dict[str, int]:
    """
    Apply offset correction to an NMR-STAR file using a LACS JSON output.
    (Applies only CA/CB/C/N offsets; H is not supported by the STAR writer here.)
    """
    if apply_selected_offsets_and_note is None:
        raise RuntimeError("apply_selected_offsets_and_note not available. Ensure apply_lacs_correction.py is on PYTHONPATH.")

    with open(lacs_output, "r", encoding="utf-8") as f:
        report = json.load(f)

    offsets_uc = _extract_offsets_for_list(report, list_id)
    atoms_use = _normalize_atoms(atoms)

    parts = [f"{a}={offsets_uc[a]:+g}" for a in atoms_use]
    details = (
        f"LACS correction applied to list_id {list_id} ({data_id}): "
        + (", ".join(parts) if parts else "none")
        + ". Source: LACS JSON."
    )

    counts = apply_selected_offsets_and_note(
        input_path=str_file,
        output_path=output_corrected,
        list_id=int(list_id),
        offsets=offsets_uc,
        atoms=atoms_use,
        release_author=release_author,
        release_details=details,
    )
    return counts


# --------------------------- types & light AA map -----------------------------

ResidueKey = Tuple[str, str, str, str]  # (Entity_ID, Entity_assembly_ID, Comp_index_ID, Comp_ID)

AA3_TO_1 = {
    'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLU':'E','GLN':'Q','GLY':'G',
    'HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S',
    'THR':'T','TRP':'W','TYR':'Y','VAL':'V',
    # common variants
    'HID':'H','HIE':'H','HIP':'H','HSD':'H','HSE':'H','HSP':'H',
    'CYX':'C','CSE':'C','CSO':'C','MSE':'M','SEC':'U','PYL':'O'
}

def tag_to_label(tag: tuple) -> str:
    """Convert a ResidueKey to a short label like ``56H`` (index+one-letter)."""
    try:
        idx = str(tag[2])
    except Exception:
        idx = "?"
    comp = (tag[3] or "").upper()
    aa1 = AA3_TO_1.get(comp, comp[:1] if comp else "?")
    return f"{idx}{aa1}"


# =============================================================================
# Core utilities
# =============================================================================

@dataclass
class FitResult:
    """Container for per-side fit results."""
    slope_pos: float
    intercept_pos: float
    fitted_pos: List[float]
    resid_pos: List[float]
    x_pos: List[float]
    y_pos: List[float]
    tags_pos: List[ResidueKey]

    slope_neg: float
    intercept_neg: float
    fitted_neg: List[float]
    resid_neg: List[float]
    x_neg: List[float]
    y_neg: List[float]
    tags_neg: List[ResidueKey]


def read_star(file_name: str) -> Dict[str, Dict[ResidueKey, Dict[str, float]]]:
    """Read an NMR-STAR file and extract atom chemical shifts."""
    try:
        ent = pynmrstar.Entry.from_file(file_name)
    except FileNotFoundError:
        return {}

    loops = ent.get_loops_by_category('Atom_chem_shift')
    if not loops:
        return {}
    col = loops[0].get_tag_names()

    idx_entity = col.index('_Atom_chem_shift.Entity_ID')
    idx_assy = col.index('_Atom_chem_shift.Entity_assembly_ID')
    idx_comp_index = col.index('_Atom_chem_shift.Comp_index_ID')
    idx_comp = col.index('_Atom_chem_shift.Comp_ID')
    idx_atom = col.index('_Atom_chem_shift.Atom_ID')
    idx_val = col.index('_Atom_chem_shift.Val')
    idx_list = col.index('_Atom_chem_shift.Assigned_chem_shift_list_ID')

    cs_data: Dict[str, Dict[ResidueKey, Dict[str, float]]] = {}
    for cs_loop in loops:
        if not cs_loop.data:
            continue
        list_id = cs_loop.data[0][idx_list]
        cs_data.setdefault(list_id, {})
        for row in cs_loop.data:
            key: ResidueKey = (row[idx_entity], row[idx_assy], row[idx_comp_index], row[idx_comp])
            cs_data[list_id].setdefault(key, {})
            try:
                cs_data[list_id][key][row[idx_atom]] = float(row[idx_val])
            except (TypeError, ValueError):
                continue
    return cs_data


def _get_random_coil(rc: RandomCoil, comp_id: str, atom: str, rc_model) -> Optional[float]:
    """Robust RC getter that tries common synonyms for amide proton."""
    try:
        return rc.get_value(comp_id, atom, rc_model)
    except Exception:
        # for H: try synonyms HN/H1
        if atom == 'H':
            for alt in ('HN', 'H1'):
                try:
                    return rc.get_value(comp_id, alt, rc_model)
                except Exception:
                    continue
        return None


def compute_deltas(resmap: Dict[ResidueKey, Dict[str, float]],
                   rc_model: Optional[Sequence[str] | str]) -> Tuple[Dict[str, List[float]], Dict[str, List[ResidueKey]]]:
    """
    Compute Δδ for H (amide) and N, and x = ΔδN − ΔδH, assembling arrays.

    Returns
    -------
    d, tags
      d contains: {'h','n','x_for_hn'}
      tags contains residue tags for 'h' and 'n'
    """
    rc = RandomCoil()
    d: Dict[str, List[float]] = {'h': [], 'n': [], 'x_for_hn': []}
    tags: Dict[str, List[ResidueKey]] = {'h': [], 'n': []}

    for residue, atom_map in resmap.items():
        comp_id = residue[-1]
        if comp_id is None or not isinstance(comp_id, str):
            continue
        try:
            # observed amide H: prefer 'H', then 'HN', then 'H1'
            obs_h = None
            for a in ('H', 'HN', 'H1'):
                if a in atom_map:
                    obs_h = float(atom_map.get(a))
                    break
            obs_n = atom_map.get('N')
            if obs_h is None or obs_n is None:
                continue

            rc_h = _get_random_coil(rc, comp_id, 'H', rc_model)
            rc_n = _get_random_coil(rc, comp_id, 'N', rc_model)
            if rc_h is None or rc_n is None:
                continue

            dd_h = obs_h - rc_h
            dd_n = float(obs_n) - rc_n
            x = dd_n - dd_h

            d['h'].append(dd_h); tags['h'].append(residue)
            d['n'].append(dd_n); tags['n'].append(residue)
            d['x_for_hn'].append(x)
        except Exception:
            continue
    return d, tags


def mad(arr: np.ndarray) -> float:
    med = np.median(arr)
    m = np.median(np.abs(arr - med))
    return float(m) if m > 0 else 1.0


def logistic_prob(z: np.ndarray, slope: float = 6.0, hinge: float = 1.0) -> np.ndarray:
    return 1.0 / (1.0 + np.exp(-slope * (z - hinge)))


def outlier_stats(residuals: np.ndarray, scale: Optional[float] = None, cutoff_k: float = 5.0) -> Tuple[List[int], List[float]]:
    if residuals.size == 0:
        return [], []
    if scale is None:
        scale = mad(residuals)
    z = np.abs(residuals) / (cutoff_k * (scale if scale > 0 else 1.0))
    flags = (z > 1.0).astype(int).tolist()
    probs = logistic_prob(z).tolist()
    return flags, probs


def _ci_and_sd(samples: np.ndarray, level: float = 0.95) -> Tuple[float, float, float]:
    if samples.size == 0:
        return float('nan'), float('nan'), float('nan')
    alpha = (1 - level) / 2.0
    lo, hi = np.quantile(samples, [alpha, 1 - alpha])
    sd = float(np.std(samples, ddof=1))
    return float(lo), float(hi), sd


def collect_and_report(fits: Dict[str, FitResult], cutoff_k: float = 5.0) -> Dict[str, Dict]:
    """Assemble offsets and outlier lists into a serializable report."""
    offsets: Dict[str, float] = {}
    offsets_split: Dict[str, Dict[str, Optional[float]]] = {}
    used_side: Dict[str, str] = {}

    for atom, fr in fits.items():
        has_pos = len(fr.x_pos) > 0
        has_neg = len(fr.x_neg) > 0
        offsets_split[atom] = {
            'pos': None if not has_pos else round(-fr.intercept_pos, 3),
            'neg': None if not has_neg else round(-fr.intercept_neg, 3),
        }
        if has_pos and has_neg:
            off = -0.5 * (fr.intercept_pos + fr.intercept_neg)
            used_side[atom] = 'both'
        elif has_pos:
            off = -fr.intercept_pos
            used_side[atom] = 'pos'
        elif has_neg:
            off = -fr.intercept_neg
            used_side[atom] = 'neg'
        else:
            off = 0.0
            used_side[atom] = 'none'
        offsets[atom] = round(float(off), 3)

    report: Dict[str, Dict] = {"offsets": offsets, "offsets_split": offsets_split, "outliers": {}, "meta": {"cutoff_k": cutoff_k, "used_side": used_side}}

    # Outliers from residuals on whichever sides were fit
    for atom, fr in fits.items():
        resid_all = np.array(fr.resid_pos + fr.resid_neg, dtype=float)
        tags_all  = fr.tags_pos + fr.tags_neg
        flags, probs = outlier_stats(resid_all, cutoff_k=cutoff_k)
        out: List[Dict[str, object]] = []
        for tag, f, p, r in zip(tags_all, flags, probs, resid_all.tolist()):
            out.append({
                "residue_key": {"entity": tag[0], "assembly": tag[1], "index": tag[2], "comp": tag[3]},
                "residual": round(float(r), 4),
                "flag": int(f),
                "prob": round(float(p), 4),
            })
        report["outliers"][atom] = out
    return report


def collect_and_report_bayes(fits: Dict[str, FitResult],
                             alpha_samples: Dict[str, Dict[str, np.ndarray]],
                             cutoff_k: float = 5.0) -> Dict[str, Dict]:
    """Assemble report including Bayesian **offset uncertainties**."""
    base = collect_and_report(fits, cutoff_k=cutoff_k)
    offsets_bayes = {}
    offsets_bayes_split = {}

    for atom, fr in fits.items():
        pos = -alpha_samples.get(atom, {}).get('pos', np.array([], float))
        neg = -alpha_samples.get(atom, {}).get('neg', np.array([], float))
        if pos.size > 0 and neg.size > 0:
            n = min(pos.size, neg.size)
            offs_draws = 0.5 * (pos[:n] + neg[:n])
        elif pos.size > 0:
            offs_draws = pos
        elif neg.size > 0:
            offs_draws = neg
        else:
            offs_draws = np.array([], float)

        mean = float(np.mean(offs_draws)) if offs_draws.size else float((fr.intercept_pos + fr.intercept_neg) / 2.0)
        lo, hi, sd = _ci_and_sd(offs_draws, level=0.95)
        offsets_bayes[atom] = {
            "mean": round(mean, 4),
            "ci95": [None if np.isnan(lo) else round(lo, 4),
                     None if np.isnan(hi) else round(hi, 4)],
            "sd": None if np.isnan(sd) else round(sd, 4),
        }

        # split stats (optional diagnostic)
        lo_p, hi_p, sd_p = _ci_and_sd(pos, level=0.95)
        offsets_bayes_split.setdefault(atom, {})['pos'] = {
            "mean": None if pos.size == 0 else round(float(np.mean(pos)), 4),
            "ci95": [None if np.isnan(lo_p) else round(lo_p, 4),
                     None if np.isnan(hi_p) else round(hi_p, 4)],
            "sd": None if np.isnan(sd_p) else round(sd_p, 4),
        }
        lo_n, hi_n, sd_n = _ci_and_sd(neg, level=0.95)
        offsets_bayes_split.setdefault(atom, {})['neg'] = {
            "mean": None if neg.size == 0 else round(float(np.mean(neg)), 4),
            "ci95": [None if np.isnan(lo_n) else round(lo_n, 4),
                     None if np.isnan(hi_n) else round(hi_n, 4)],
            "sd": None if np.isnan(sd_n) else round(sd_n, 4),
        }

    base["offsets_bayes"] = offsets_bayes
    base["offsets_bayes_split"] = offsets_bayes_split
    return base


# =============================================================================
# Plotting (optional)
# =============================================================================

def _plot_atom(atom: str, fr: FitResult, outdir: Path, data_id: str, list_id, method: str, cutoff_k: float = 5.0) -> None:
    """Render interactive plots for a single nucleus with residue labels."""
    if not PLOTLY_AVAILABLE:
        return

    x_pos, y_pos = np.asarray(fr.x_pos, float), np.asarray(fr.y_pos, float)
    x_neg, y_neg = np.asarray(fr.x_neg, float), np.asarray(fr.y_neg, float)

    labels_pos = [tag_to_label(t) for t in fr.tags_pos]
    labels_neg = [tag_to_label(t) for t in fr.tags_neg]

    def line_xy(x, slope, intercept):
        if x.size == 0:
            return np.array([]), np.array([])
        xs = np.asarray(sorted(x.tolist()), float)
        ys = intercept + slope * xs
        return xs, ys

    xlp, ylp = line_xy(x_pos, fr.slope_pos, fr.intercept_pos)
    xln, yln = line_xy(x_neg, fr.slope_neg, fr.intercept_neg)

    resid_all = np.array(fr.resid_pos + fr.resid_neg, dtype=float)
    flags, probs = outlier_stats(resid_all, cutoff_k=cutoff_k)
    n_pos = len(fr.resid_pos)
    flags_pos, flags_neg = flags[:n_pos], flags[n_pos:]
    probs_pos, probs_neg = probs[:n_pos], probs[n_pos:]

    fig = go.Figure()

    def add_side(x, y, labels, flags, probs, name_in, name_out, line_x, line_y, dashed=False):
        if x.size:
            flags_arr = np.asarray(flags, int) if flags else np.zeros(x.shape[0], dtype=int)
            probs_arr = np.asarray(probs, float) if probs else np.zeros(x.shape[0], dtype=float)
            mask_out = flags_arr == 1
            mask_in  = ~mask_out

            if mask_in.any():
                fig.add_trace(go.Scatter(
                    x=x[mask_in], y=y[mask_in], mode="markers+text",
                    name=name_in,
                    text=[labels[i] for i in np.where(mask_in)[0]],
                    textposition="top center",
                    textfont=dict(size=10),
                    customdata=np.c_[flags_arr[mask_in], probs_arr[mask_in], np.where(mask_in)[0]],
                    hovertemplate=(
                        "x=%{x:.3f}<br>y=%{y:.3f}<br>label=%{text}"
                        "<br>outlier=%{customdata[0]}<br>prob=%{customdata[1]:.2f}"
                    ),
                    marker=dict(symbol="circle", size=8),
                ))
            if mask_out.any():
                fig.add_trace(go.Scatter(
                    x=x[mask_out], y=y[mask_out], mode="markers+text",
                    name=name_out,
                    text=[labels[i] for i in np.where(mask_out)[0]],
                    textposition="top center",
                    textfont=dict(size=10),
                    customdata=np.c_[flags_arr[mask_out], probs_arr[mask_out], np.where(mask_out)[0]],
                    hovertemplate=(
                        "x=%{x:.3f}<br>y=%{y:.3f}<br>label=%{text}"
                        "<br>outlier=%{customdata[0]}<br>prob=%{customdata[1]:.2f}"
                    ),
                    marker=dict(symbol="x", size=12, line=dict(width=1.5)),
                ))
            if line_x.size:
                fig.add_trace(go.Scatter(
                    x=line_x, y=line_y, mode="lines",
                    name=f"{'fit +' if not dashed else 'fit −'}",
                    line=dict(dash="dash" if dashed else "solid")
                ))

    add_side(x_pos, y_pos, labels_pos, flags_pos, probs_pos, "x ≥ 0 inliers", "x ≥ 0 outliers", xlp, ylp, dashed=False)
    add_side(x_neg, y_neg, labels_neg, flags_neg, probs_neg, "x < 0 inliers", "x < 0 outliers", xln, yln, dashed=True)

    ytitle = {'h':'ΔδH','n':'ΔδN'}.get(atom, 'Δδ')
    x_label = "ΔδN − ΔδH"
    fig.update_layout(
        title=f"{data_id}: {atom.upper()} vs {x_label}",
        xaxis_title=x_label,
        yaxis_title=ytitle,
        legend_orientation="h",
        template="plotly_white"
    )
    outdir.mkdir(parents=True, exist_ok=True)
    html_path = outdir / f"{data_id}_{atom}_{list_id}_{method}.html"
    fig.write_html(html_path, include_mathjax="cdn")

    # Probability bar
    def safe_res_index(tag):
        try:
            return int(tag[2])
        except Exception:
            return None

    res_indices = [safe_res_index(t) for t in (fr.tags_pos + fr.tags_neg)]
    xvals = [ri if ri is not None else i+1 for i, ri in enumerate(res_indices)]
    fig2 = px.bar(x=xvals, y=probs, labels={"x": "Residue Index", "y": "Outlier probability"})
    fig2.update_layout(title=f"{data_id}: {atom.upper()} outlier probability", yaxis_range=[0,1], template="plotly_white")
    fig2.write_html(outdir / f"{data_id}_{atom}_{list_id}_{method}_prob.html")

    try:  # optional PDF
        fig.write_image(outdir / f"{data_id}_{atom}_{list_id}_{method}.pdf")
        fig2.write_image(outdir / f"{data_id}_{atom}_{list_id}_{method}prob.pdf")
    except Exception:
        pass


def maybe_plot_all(fits: Dict[str, FitResult], outdir: Optional[Path], data_id: str, method: str, enable_plots: bool, list_id, cutoff_k: float = 5.0) -> None:
    if not enable_plots or not fits:
        return
    if not PLOTLY_AVAILABLE:  # pragma: no cover - optional
        print("[plot] plotly not available; skipping plots.", file=sys.stderr)
        return
    out = Path(outdir or (Path.cwd() / "lacs_output"))
    out.mkdir(parents=True, exist_ok=True)
    for atom, fr in fits.items():
        _plot_atom(atom, fr, out, data_id, list_id, method=method, cutoff_k=cutoff_k)


# =============================================================================
# Method-specific fitters (unchanged)
# =============================================================================

def fit_side_rlm_tukey(x: Sequence[float], y: Sequence[float]) -> Tuple[float, float, List[float], List[float]]:
    try:
        import statsmodels.api as sm
        from statsmodels.robust.norms import TukeyBiweight
    except Exception as e:
        raise SystemExit("This method requires 'statsmodels' (pip install statsmodels).") from e

    X = sm.add_constant(np.asarray(x, dtype=float))
    yv = np.asarray(y, dtype=float)
    model = sm.RLM(yv, X, M=TukeyBiweight(c=4.685))
    try:
        res = model.fit(scale_est="mad")
    except Exception:
        res = model.fit()

    fitted = res.fittedvalues.astype(float).tolist()
    resid = (yv - res.fittedvalues).astype(float).tolist()
    slope = float(res.params[1])
    intercept = float(res.params[0])
    return slope, intercept, fitted, resid


def fit_side_theilsen(x: Sequence[float], y: Sequence[float]) -> Tuple[float, float, List[float], List[float]]:
    try:
        from sklearn.linear_model import TheilSenRegressor
    except Exception as e:
        raise SystemExit("This method requires scikit-learn (pip install scikit-learn).") from e

    X = np.asarray(x, dtype=float).reshape(-1,1)
    yv = np.asarray(y, dtype=float)
    reg = TheilSenRegressor(random_state=0).fit(X, yv)
    yhat = reg.predict(X)
    resid = (yv - yhat).astype(float).tolist()
    slope = float(reg.coef_[0])
    intercept = float(reg.intercept_)
    return slope, intercept, yhat.astype(float).tolist(), resid


def fit_side_ransac(x: Sequence[float], y: Sequence[float]) -> Tuple[float, float, List[float], List[float]]:
    try:
        from sklearn.linear_model import RANSACRegressor, LinearRegression
    except Exception as e:
        raise SystemExit("This method requires scikit-learn (pip install scikit-learn).") from e

    X = np.asarray(x, dtype=float).reshape(-1,1)
    yv = np.asarray(y, dtype=float)

    lr0 = LinearRegression().fit(X, yv)
    resid0 = yv - lr0.predict(X)
    mad0 = np.median(np.abs(resid0 - np.median(resid0)))
    sigma = (1.4826 * mad0) if mad0 > 0 else (np.std(resid0) if resid0.size > 1 else 1.0)
    thresh = 2.5 * sigma

    ransac = RANSACRegressor(
        estimator=LinearRegression(),
        residual_threshold=thresh,
        random_state=0,
    ).fit(X, yv)

    yhat = ransac.predict(X)
    resid = (yv - yhat).astype(float).tolist()
    slope = float(ransac.estimator_.coef_[0])
    intercept = float(ransac.estimator_.intercept_)
    return slope, intercept, yhat.astype(float).tolist(), resid


def fit_side_quantile(x: Sequence[float], y: Sequence[float]) -> Tuple[float, float, List[float], List[float]]:
    try:
        import statsmodels.api as sm
    except Exception as e:
        raise SystemExit("This method requires 'statsmodels' (pip install statsmodels).") from e

    X = sm.add_constant(np.asarray(x, dtype=float))
    yv = np.asarray(y, dtype=float)
    mod = sm.QuantReg(yv, X)
    res = mod.fit(q=0.5)
    fitted_line = (res.params[0] + res.params[1] * X[:,1]).astype(float)
    resid = (yv - fitted_line).astype(float).tolist()
    slope = float(res.params[1])
    intercept = float(res.params[0])
    return slope, intercept, fitted_line.tolist(), resid


def fit_side_bayes_t(x: Sequence[float], y: Sequence[float]) -> Tuple[float, float, List[float], List[float], float, float, np.ndarray]:
    try:
        import pymc as pm
    except Exception as e:
        raise SystemExit("This method requires PyMC (pip install pymc).") from e

    x_arr = np.asarray(x, dtype=float)
    y_arr = np.asarray(y, dtype=float)
    with pm.Model() as m:
        alpha = pm.Normal("alpha", 0, 10)
        beta  = pm.Normal("beta", 0, 10)
        sigma = pm.HalfNormal("sigma", 2.0)
        nu    = pm.Exponential("nu_raw", 1/30) + 1  # df > 1
        mu = alpha + beta * x_arr
        pm.StudentT("y", nu=nu, mu=mu, sigma=sigma, observed=y_arr)
        idata = pm.sample(draws=1000, tune=1000, target_accept=0.9, chains=4, cores=4, progressbar=False, random_seed=42)
        post = idata.posterior

    def _as_values(a):
        try:
            return a.values
        except AttributeError:
            try:
                return a.to_numpy()
            except AttributeError:
                return np.asarray(a)

    alpha_arr = _as_values(post["alpha"])
    beta_arr  = _as_values(post["beta"])
    sigma_arr = _as_values(post["sigma"])
    nu_arr    = _as_values(post["nu_raw"])

    alpha_mean = float(np.mean(alpha_arr))
    beta_mean  = float(np.mean(beta_arr))
    sigma_mean = float(np.mean(sigma_arr))
    nu_mean    = float(np.mean(nu_arr) + 1.0)

    alpha_draws = alpha_arr.reshape(-1)

    fitted = (alpha_mean + beta_mean * x_arr).astype(float).tolist()
    resid = (y_arr - (alpha_mean + beta_mean * x_arr)).astype(float).tolist()
    return beta_mean, alpha_mean, fitted, resid, nu_mean, sigma_mean, alpha_draws


# ---------------------- nucleus-level orchestration ---------------------------

def _fit_atom_by_method(method: str, xvals: List[float], yvals: List[float], tags: List[ResidueKey],
                        min_per_side: int = MIN_PER_SIDE_DEFAULT) -> FitResult:
    """Fit a single nucleus by splitting x by sign and applying ``method``."""
    x = np.asarray(xvals, dtype=float); y = np.asarray(yvals, dtype=float)
    pos = x >= 0; neg = ~pos

    def do_side(mask):
        if int(mask.sum()) < int(min_per_side):
            # Not enough data for this side; mark as empty
            return 1.0, 0.0, [], [], [], [], []
        xx, yy = x[mask], y[mask]
        if method == "tukey":
            s, b, f, r = fit_side_rlm_tukey(xx, yy)
        elif method == "theilsen":
            s, b, f, r = fit_side_theilsen(xx, yy)
        elif method == "ransac":
            s, b, f, r = fit_side_ransac(xx, yy)
        elif method == "quantile":
            s, b, f, r = fit_side_quantile(xx, yy)
        elif method == "bayes":
            s, b, f, r, _, _, _ = fit_side_bayes_t(xx, yy)
        else:
            raise ValueError(f"Unknown method: {method}")
        tg = [t for t, m in zip(tags, mask) if m]
        return s, b, f, r, xx.tolist(), yy.tolist(), tg

    s_p, b_p, f_p, r_p, xp, yp, t_p = do_side(pos)
    s_n, b_n, f_n, r_n, xn, yn, t_n = do_side(neg)

    return FitResult(s_p, b_p, f_p, r_p, xp, yp, t_p, s_n, b_n, f_n, r_n, xn, yn, t_n)


def _fit_atom_bayes(xvals: List[float], yvals: List[float], tags: List[ResidueKey],
                    min_per_side: int = MIN_PER_SIDE_DEFAULT) -> Tuple[FitResult, Dict[str, np.ndarray]]:
    """Bayesian version of :func:`_fit_atom_by_method` that returns intercept draws."""
    x = np.asarray(xvals, dtype=float); y = np.asarray(yvals, dtype=float)
    pos = x >= 0; neg = ~pos
    draws: Dict[str, np.ndarray] = {"pos": np.array([], float), "neg": np.array([], float)}

    def do_side(mask, side_key):
        if int(mask.sum()) < int(min_per_side):
            return 1.0, 0.0, [], [], [], [], []
        xx, yy = x[mask], y[mask]
        s, b, f, r, _, _, alpha_draws = fit_side_bayes_t(xx, yy)
        draws[side_key] = alpha_draws
        tg = [t for t, m in zip(tags, mask) if m]
        return s, b, f, r, xx.tolist(), yy.tolist(), tg

    s_p, b_p, f_p, r_p, xp, yp, t_p = do_side(pos, "pos")
    s_n, b_n, f_n, r_n, xn, yn, t_n = do_side(neg, "neg")

    fr = FitResult(s_p, b_p, f_p, r_p, xp, yp, t_p, s_n, b_n, f_n, r_n, xn, yn, t_n)
    return fr, draws


def run_lacs(str_file: str, method: str, data_id: str, rc_model: Optional[Sequence[str] | str] = None,
             outdir: Optional[Path] = None, plots: bool = True, cutoff_k: float = 5.0,
             min_per_side: int = MIN_PER_SIDE_DEFAULT) -> Dict[str, Dict]:
    """
    Run the selected robust method over an NMR-STAR file for H and N using x=ΔδN−ΔδH.
    Returns a report keyed by chemical-shift list id.
    """
    cs = read_star(str_file)
    results: Dict[str, Dict] = {}
    for list_id, resmap in cs.items():
        d, tags = compute_deltas(resmap, rc_model)
        fits: Dict[str, FitResult] = {}
        alpha_samples: Dict[str, Dict[str, np.ndarray]] = {}

        # Only H and N; both use the same x key
        for atom in ('h', 'n'):
            yvals = d[atom]
            xvals = d['x_for_hn']
            tg = tags[atom]
            if len(yvals) >= 2 and len(xvals) == len(yvals):
                if method == "bayes":
                    fr, draws = _fit_atom_bayes(xvals, yvals, tg, min_per_side=min_per_side)
                    alpha_samples[atom] = draws
                else:
                    fr = _fit_atom_by_method(method, xvals, yvals, tg, min_per_side=min_per_side)
                fits[atom] = fr

        if fits:
            maybe_plot_all(fits, outdir, data_id, method, plots, list_id, cutoff_k=cutoff_k)
            if method == "bayes":
                results[list_id] = collect_and_report_bayes(fits, alpha_samples, cutoff_k=cutoff_k)
            else:
                results[list_id] = collect_and_report(fits, cutoff_k=cutoff_k)
    return results


# =============================================================================
# CLI
# =============================================================================

def _default_json_path(outdir: Optional[Path], data_id: str, method: str) -> Path:
    base = Path.cwd() if outdir is None else Path(outdir)
    return base / f"{data_id}_{method}.json"


def main(argv=None) -> None:
    import argparse

    p = argparse.ArgumentParser(description="LACS H/N variant: robust fits for ΔδH, ΔδN vs (ΔδN−ΔδH).")
    p.add_argument("star_file", help="Path to NMR-STAR .str file")
    p.add_argument("--method", default='bayes',
                   choices=["tukey","theilsen","ransac","quantile","bayes"],
                   help="Robust regression method to use (default Bayes)")
    p.add_argument("--data-id", default="LACS_HN")
    p.add_argument("--rc-model", nargs="*", default=None,
                   help="Random-coil model alias(es), e.g. wis wan; omit for average of all")
    p.add_argument("--out", type=Path, default=None, help="Output directory for plots")
    p.add_argument("--no-plots", action="store_true", help="Disable plotting entirely")
    p.add_argument('--min-per-side', type=int, default=MIN_PER_SIDE_DEFAULT,
                   help='Minimum number of points required on each sign side (default: 5).')
    p.add_argument("--cutoff-k", type=float, default=5.0, help="Outlier cutoff multiplier (default 5.0)")
    p.add_argument("--json-out", type=Path, default=None, help="Where to write JSON (defaults to <data_id>_<method>.json)")

    # Offset application kept for CA/CB/C/N; useful if you still want to correct those from another run
    p.add_argument("--apply-offsets", action="store_true",
                   help="After computing offsets, apply them to the input .str (CA/CB/C/N only).")
    p.add_argument("--output-corrected", type=Path,
                   help="Path to write the corrected NMR-STAR file (required with --apply-offsets).")
    p.add_argument("--atoms", nargs="+",
                   choices=["CA", "CB", "C", "N", "ca", "cb", "c", "n"],
                   default=["CA", "CB", "C", "N"],
                   help="Subset of atoms to correct (default: CA CB C N).")
    p.add_argument("--release-author", default="BMRB",
                   help="Author to record in the Release loop (default: BMRB).")

    args = p.parse_args(argv)
    rc_model = None if args.rc_model == [] else (args.rc_model if args.rc_model is not None else None)

    report = run_lacs(args.star_file, method=args.method, data_id=args.data_id,
                      rc_model=rc_model, outdir=args.out, plots=not args.no_plots,
                      cutoff_k=args.cutoff_k, min_per_side=args.min_per_side)

    json_path = args.json_out or _default_json_path(args.out, args.data_id, args.method)
    json_path.parent.mkdir(parents=True, exist_ok=True)
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2)

    if args.apply_offsets:
        if apply_selected_offsets_and_note is None:
            raise SystemExit("Cannot apply offsets: apply_lacs_correction.py not importable.")
        if not args.output_corrected:
            raise SystemExit("--output-corrected is required with --apply-offsets.")

        for list_id in report:
            try:
                offsets_uc = _extract_offsets_for_list(report, list_id)
            except KeyError:
                with open(json_path, "r", encoding="utf-8") as f:
                    report2 = json.load(f)
                offsets_uc = _extract_offsets_for_list(report2, list_id)

            atoms_use = _normalize_atoms(args.atoms)
            parts = [f"{a}={offsets_uc[a]:+g}" for a in atoms_use]
            details = (
                f"LACS correction applied to list_id {list_id} ({args.data_id}): "
                + (", ".join(parts) if parts else "none")
                + ". Source: LACS (same run)."
            )
            # Write first list from input → output, subsequent lists in-place to output
            list_ids = list(report.keys())
            first = (list_ids.index(list_id) == 0)
            counts = apply_selected_offsets_and_note(
                input_path=(args.star_file if first else str(args.output_corrected)),
                output_path=str(args.output_corrected),
                list_id=int(list_id),
                offsets=offsets_uc,
                atoms=atoms_use,
                release_author=args.release_author,
                release_details=details,
            )

        print("Applied offsets to file:")
        for a in atoms_use:
            print(f"  {a}: {counts[a]}")
        print(f"Total updated: {counts['total']}")
        print(f"Wrote corrected file: {args.output_corrected}")


if __name__ == "__main__":
    main()
