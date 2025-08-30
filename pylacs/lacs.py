
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
LACS unified CLI
================

A single command-line tool to compute LACS offsets, outliers, and soft
probabilities from an NMR-STAR ``.str`` file using one of five robust
linear-fit methods:

- ``tukey``: Tukey biweight RLM (statsmodels)
- ``theilsen``: Theil–Sen median-of-slopes (scikit-learn)
- ``ransac``: RANSAC with robust threshold (scikit-learn)
- ``quantile``: Quantile regression at τ=0.5 (statsmodels)
- ``bayes``: Bayesian Student-t regression (PyMC)

Each method fits per nucleus (ΔδC, ΔδN, ΔδH, ΔδCA, ΔδCB) as

    y = b + m·x,  where x = ΔδCA − ΔδCB,

splitting the data by the sign of ``x`` (x ≥ 0 vs x < 0) and fitting a line
to each side. Offsets are reported as the average intercept across sides.
Outliers are flagged from robust-scaled residuals with a smooth probability map.

Bayesian method uncertainty
---------------------------
When ``--method bayes`` is used, the JSON additionally contains **offset
uncertainties** computed from the posterior intercept samples on each side.
We form samples of

    offset = 0.5 * (b_pos + b_neg)

and report the mean, 95% credible interval, and standard deviation per nucleus.

Usage
-----
.. code-block:: bash

    python lacs_unified.py ENTRY.str --method tukey --data-id myprot --out figs --json-out myprot_tukey.json

Dependencies
------------
- *Common:* ``pynmrstar``, your ``random_coil_refactored.py`` on PYTHONPATH
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
    from .random_coil import RandomCoil
except Exception as e:
    try:
        from random_coil import RandomCoil
    except Exception as e:
        raise SystemExit("Couldn't import random_coil.py. Ensure it is on PYTHONPATH or in the same folder.") from e

ResidueKey = Tuple[str, str, str, str]  # (Entity_ID, Entity_assembly_ID, Comp_index_ID, Comp_ID)

# Map 3-letter residue codes to one-letter; includes common variants
AA3_TO_1 = {
    'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLU':'E','GLN':'Q','GLY':'G',
    'HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S',
    'THR':'T','TRP':'W','TYR':'Y','VAL':'V',
    # common alt names / modified residues mapping to canonical letters
    'HID':'H','HIE':'H','HIP':'H','HSD':'H','HSE':'H','HSP':'H',
    'CYX':'C','CSE':'C','CSO':'C','MSE':'M','SEC':'U','PYL':'O'
}
CA_SCALE = {'ILE': 4.8, 'GLN': 4.5, 'GLY': 1.7, 'GLU': 4.5, 'CYS': 6.6, 'ASP': 4.2, 'SER': 4.1, 'LYS': 4.6, 'PRO': 3.1, 'ASN': 3.6, 'VAL': 5.8, 'THR': 4.9, 'HIS': 4.0, 'TRP': 4.0, 'PHE': 5.3, 'ALA': 5.1, 'MET': 4.2, 'LEU': 4.3, 'ARG': 4.6, 'TYR': 4.7}
CB_SCALE = {'ILE': 2.4, 'GLN': 3.5, 'GLY': 1.0, 'GLU': 3.4, 'CYS': 3.6, 'ASP': 1.9, 'SER': 2.55, 'LYS': 2.4, 'PRO': 1.0, 'ASN': 1.7, 'VAL': 2.7, 'THR': 1.5, 'HIS': 3.4, 'TRP': 1.9, 'PHE': 2.5, 'ALA': 4.2, 'MET': 2.8, 'LEU': 2.65, 'ARG': 2.7, 'TYR': 3.1}
POS_SCALE = {'ILE': 4.83, 'GLN': 3.82, 'GLY': 1.0, 'GLU': 3.24, 'CYS': 4.35, 'ASP': 3.42, 'SER': 3.98, 'LYS': 3.42, 'PRO': 3.12, 'ASN': 2.94, 'VAL': 5.19, 'THR': 5.56, 'HIS': 4.26, 'TRP': 3.09, 'PHE': 4.04, 'ALA': 4.01, 'MET': 3.44, 'LEU': 3.47, 'ARG': 3.65, 'TYR': 3.84}
NEG_SCALE = {'ILE': 2.37, 'GLN': 4.18, 'GLY': 1.0, 'GLU': 4.66, 'CYS': 5.85, 'ASP': 2.68, 'SER': 2.67, 'LYS': 3.58, 'PRO': 0.98, 'ASN': 2.36, 'VAL': 3.31, 'THR': 0.84, 'HIS': 3.14, 'TRP': 2.81, 'PHE': 3.76, 'ALA': 5.29, 'MET': 3.56, 'LEU': 3.48, 'ARG': 3.65, 'TYR': 3.96}

def tag_to_label(tag: tuple) -> str:
    """Convert a ResidueKey to a short label like ``56H`` (index+one-letter).

    - Uses ``Comp_index_ID`` as the residue number
    - Converts ``Comp_ID`` (3-letter) to a one-letter code when possible
    """
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
    """Container for per-side fit results.

    Parameters
    ----------
    slope_pos, intercept_pos : float
        Slope and intercept for the ``x ≥ 0`` side.
    fitted_pos, resid_pos : list[float]
        Fitted values and residuals for points on the positive side.
    x_pos, y_pos : list[float]
        Observed x and y values on the positive side.
    tags_pos : list[ResidueKey]
        Residue identifiers for positive-side points.

    slope_neg, intercept_neg : float
        Slope and intercept for the ``x < 0`` side.
    fitted_neg, resid_neg : list[float]
        Fitted values and residuals for points on the negative side.
    x_neg, y_neg : list[float]
        Observed x and y values on the negative side.
    tags_neg : list[ResidueKey]
        Residue identifiers for negative-side points.
    """
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
    """Read an NMR-STAR file and extract atom chemical shifts.

    Parameters
    ----------
    file_name : str
        Path to a ``.str`` file.

    Returns
    -------
    dict
        Mapping ``{list_id -> {(entity, assembly, comp_index, comp_id) -> {atom -> value}}}``.

    Notes
    -----
    Only rows in the ``Atom_chem_shift`` category are parsed. Non-numeric
    values are skipped.
    """
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


def compute_deltas(resmap: Dict[ResidueKey, Dict[str, float]], rc_model: Optional[Sequence[str] | str]) -> Tuple[Dict[str, List[float]], Dict[str, List[ResidueKey]]]:
    """Compute Δδ against random-coil and assemble arrays for each nucleus.

    Parameters
    ----------
    resmap : dict
        Chemical-shift map for a single list id as returned by :func:`read_star`.
    rc_model : sequence[str] | str | None
        Random-coil model alias(es) to pass to :class:`RandomCoil`. If ``None``,
        the model averages all available references.

    Returns
    -------
    (dict, dict)
        The first dict holds numeric arrays:
        ``{'ca','cb','c','n','ha','x_for_c','x_for_n','x_for_h'}``.
        The second dict holds residue tags per nucleus.

    Notes
    -----
    - ``x = ΔδCA − ΔδCB`` is computed only when both CA and CB are present.
    - Δδ values are formed by subtracting the random-coil reference for the
      residue type from the observed chemical shift.
    """
    rc = RandomCoil()
    d: Dict[str, List[float]] = {
        'ca': [], 'cb': [], 'c': [], 'n': [], 'ha': [],
        'x_for_c': [], 'x_for_n': [], 'x_for_ha': [],
    }
    tags: Dict[str, List[ResidueKey]] = {k: [] for k in ['ca','cb','c','n','ha']}

    for residue, atom_map in resmap.items():
        comp_id = residue[-1]
        try:
            wca = CA_SCALE[comp_id]
        except KeyError:
            wca = 1.0
        try:
            wcb = CB_SCALE[comp_id]
        except KeyError:
            wcb = 1.0
        if comp_id is None or not isinstance(comp_id, str):
            continue
        try:
            ca = atom_map.get('CA')
            cb = atom_map.get('CB')
            c  = atom_map.get('C')
            n  = atom_map.get('N')
            ha  = atom_map.get('HA')

            if ca is not None: ca = float(ca) - rc.get_value(comp_id, 'CA', rc_model)
            if cb is not None: cb = float(cb) - rc.get_value(comp_id, 'CB', rc_model)
            if c  is not None: c  = float(c)  - rc.get_value(comp_id, 'C',  rc_model)
            if n  is not None: n  = float(n)  - rc.get_value(comp_id, 'N',  rc_model)
            if ha  is not None: ha  = float(ha)  - rc.get_value(comp_id, 'HA',  rc_model)

            if ca is not None and cb is not None:
                x = ca - cb
                # if x1 >=0 :
                #     x=x1/POS_SCALE[comp_id]
                # else:
                #     x=x1/NEG_SCALE[comp_id]


                d['ca'].append(ca); tags['ca'].append(residue)
                d['cb'].append(cb); tags['cb'].append(residue)
                if c is not None:
                    d['c'].append(c); d['x_for_c'].append(x); tags['c'].append(residue)
                if n is not None:
                    d['n'].append(n); d['x_for_n'].append(x); tags['n'].append(residue)
                if ha is not None:
                    d['ha'].append(ha); d['x_for_ha'].append(x); tags['ha'].append(residue)
        except Exception:
            continue
    return d, tags


def mad(arr: np.ndarray) -> float:
    """Median absolute deviation (MAD).

    Parameters
    ----------
    arr : ndarray

    Returns
    -------
    float
        MAD of the input, or 1.0 if the MAD is zero (prevents divide-by-zero).
    """
    med = np.median(arr)
    m = np.median(np.abs(arr - med))
    return float(m) if m > 0 else 1.0


def logistic_prob(z: np.ndarray, slope: float = 6.0, hinge: float = 1.0) -> np.ndarray:
    """Logistic map from a nonnegative score to a soft probability in [0, 1].

    Parameters
    ----------
    z : ndarray
        Nonnegative score (e.g., robust z-score), where 1 corresponds to the
        desired cutoff.
    slope : float, optional
        Steepness of the logistic curve.
    hinge : float, optional
        Center point of the curve (probability 0.5 at ``z = hinge``).

    Returns
    -------
    ndarray
        Probabilities of the same shape as ``z``.
    """
    return 1.0 / (1.0 + np.exp(-slope * (z - hinge)))


def outlier_stats(residuals: np.ndarray, scale: Optional[float] = None, cutoff_k: float = 5.0) -> Tuple[List[int], List[float]]:
    """Compute 0/1 outlier flags and smooth probabilities from residuals.

    Parameters
    ----------
    residuals : ndarray
        Fitted residuals.
    scale : float, optional
        Robust scale (MAD) to use. If ``None``, MAD is computed from data.
    cutoff_k : float, optional
        Scale multiplier for the cutoff (default 5.0). Points with
        ``|r|/(k·MAD) > 1`` are flagged as outliers.

    Returns
    -------
    (list[int], list[float])
        Flags and probabilities for each residual.
    """
    if residuals.size == 0:
        return [], []
    if scale is None:
        scale = mad(residuals)
    z = np.abs(residuals) / (cutoff_k * (scale if scale > 0 else 1.0))
    flags = (z > 1.0).astype(int).tolist()
    probs = logistic_prob(z).tolist()
    return flags, probs


def _ci_and_sd(samples: np.ndarray, level: float = 0.95) -> Tuple[float, float, float]:
    """Compute central credible interval and standard deviation.

    Parameters
    ----------
    samples : ndarray
        1D array of posterior draws.
    level : float
        Credible interval mass (default 0.95).

    Returns
    -------
    (lo, hi, sd) : tuple[float, float, float]
        Lower and upper bounds of the central credible interval and the SD.
    """
    if samples.size == 0:
        return float('nan'), float('nan'), float('nan')
    alpha = (1 - level) / 2.0
    lo, hi = np.quantile(samples, [alpha, 1 - alpha])
    sd = float(np.std(samples, ddof=1))
    return float(lo), float(hi), sd


def collect_and_report(fits: Dict[str, FitResult], cutoff_k: float = 5.0) -> Dict[str, Dict]:
    """Assemble offsets and outlier lists into a serializable report.

    Parameters
    ----------
    fits : dict[str, FitResult]
        Per-nucleus fit results.
    cutoff_k : float
        Cutoff multiplier used for outlier scoring (for reproducibility in logs).

    Returns
    -------
    dict
        Report with keys ``offsets`` and ``outliers``.
    """
    offsets = {atom: round((fr.intercept_pos + fr.intercept_neg) / 2.0, 3) for atom, fr in fits.items()}
    report: Dict[str, Dict] = {"offsets": offsets, "outliers": {}, "meta": {"cutoff_k": cutoff_k}}
    for atom, fr in fits.items():
        resid_all = np.array(fr.resid_pos + fr.resid_neg, dtype=float)
        tags_all  = fr.tags_pos + fr.tags_neg
        flags, probs = outlier_stats(resid_all, cutoff_k=cutoff_k)
        out = []
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
                             cutoff_k: float = 2.5) -> Dict[str, Dict]:
    """Assemble report including Bayesian **offset uncertainties**.

    Parameters
    ----------
    fits : dict[str, FitResult]
        Per-nucleus fit results (means).
    alpha_samples : dict[str, dict[str, ndarray]]
        Posterior draws of per-side intercepts for each nucleus.
        Structure: ``alpha_samples[atom]['pos']`` and ``['neg']`` -> 1D arrays.
    cutoff_k : float
        Outlier cutoff multiplier used for residual-based probabilities.

    Returns
    -------
    dict
        Report with keys:
        - ``offsets``: mean offsets (as before)
        - ``offsets_bayes``: dict with ``mean``, ``ci95``, ``sd`` per atom
        - ``outliers``: list of outlier diagnostics per atom
        - ``meta.cutoff_k``: scalar
    """
    base = collect_and_report(fits, cutoff_k=cutoff_k)
    offsets_bayes = {}
    for atom, fr in fits.items():
        pos = alpha_samples.get(atom, {}).get('pos', np.array([], float))
        neg = alpha_samples.get(atom, {}).get('neg', np.array([], float))
        if pos.size > 0 and neg.size > 0:
            n = min(pos.size, neg.size)
            offs_draws = 0.5 * (pos[:n] + neg[:n])
        elif pos.size > 0:
            offs_draws = pos
        elif neg.size > 0:
            offs_draws = neg
        else:
            offs_draws = np.array([], float)
        mean = float(np.mean(offs_draws)) if offs_draws.size else float((fr.intercept_pos + fr.intercept_neg)/2.0)
        lo, hi, sd = _ci_and_sd(offs_draws, level=0.95)
        offsets_bayes[atom] = {
            "mean": round(mean, 4),
            "ci95": [None if np.isnan(lo) else round(lo, 4),
                     None if np.isnan(hi) else round(hi, 4)],
            "sd": None if np.isnan(sd) else round(sd, 4),
        }
    base["offsets_bayes"] = offsets_bayes
    return base


# =============================================================================
# Plotting (optional): scatter with outlier marks and probability bars
# =============================================================================


def _plot_atom(atom: str, fr: FitResult, outdir: Path, data_id: str, method: str, cutoff_k: float = 5.0) -> None:
    """Render interactive plots for a single nucleus with residue labels.

    Parameters
    ----------
    atom : str
        Nucleus key (``'ca','cb','c','n','ha'``).
    fr : FitResult
        Fit results for the nucleus.
    outdir : Path
        Output directory for plots (HTML and optional PDF via kaleido).
    data_id : str
        Identifier to include in figure titles/filenames.
    cutoff_k : float, optional
        Outlier cutoff multiplier used when coloring points.
    """
    if not PLOTLY_AVAILABLE:
        return

    x_pos, y_pos = np.asarray(fr.x_pos, float), np.asarray(fr.y_pos, float)
    x_neg, y_neg = np.asarray(fr.x_neg, float), np.asarray(fr.y_neg, float)

    # Pre-compute labels like "56H"
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

    def add_side(x, y, labels, tags, flags, probs, name_in, name_out, line_x, line_y, dashed=False):
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

    add_side(x_pos, y_pos, labels_pos, fr.tags_pos, flags_pos, probs_pos, "x ≥ 0 inliers", "x ≥ 0 outliers", xlp, ylp, dashed=False)
    add_side(x_neg, y_neg, labels_neg, fr.tags_neg, flags_neg, probs_neg, "x < 0 inliers", "x < 0 outliers", xln, yln, dashed=True)

    ytitle = {'ca':'ΔδCA','cb':'ΔδCB','c':'ΔδC','n':'ΔδN','ha':'ΔδHA'}.get(atom, 'Δδ')
    fig.update_layout(
        title=f"{data_id}: {atom.upper()} vs ΔδCA−ΔδCB",
        xaxis_title="ΔδCA − ΔδCB",
        yaxis_title=ytitle,
        legend_orientation="h",
        template="plotly_white"
    )
    outdir.mkdir(parents=True, exist_ok=True)
    html_path = outdir / f"{data_id}_{atom}_{method}.html"
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
    fig2.write_html(outdir / f"{data_id}_{atom}_{method}_prob.html")

    try:  # optional PDF
        fig.write_image(outdir / f"{data_id}_{atom}_{method}.pdf")
        fig2.write_image(outdir / f"{data_id}_{atom}_{method}prob.pdf")
    except Exception:
        pass


def maybe_plot_all(fits: Dict[str, FitResult], outdir: Optional[Path], data_id: str, method: str, enable_plots: bool, cutoff_k: float = 5.0) -> None:
    """Render all per-nucleus plots if plotting is enabled and supported.

    Parameters
    ----------
    fits : dict[str, FitResult]
        Per-nucleus fit results.
    outdir : Path | None
        Directory for plot outputs. If ``None``, ``./lacs_output`` is used.
    data_id : str
        Identifier included in figure titles/filenames.
    enable_plots : bool
        If ``False`` or Plotly is unavailable, plotting is skipped.
    cutoff_k : float, optional
        Outlier cutoff multiplier used when coloring points.
    """
    if not enable_plots or not fits:
        return
    if not PLOTLY_AVAILABLE:  # pragma: no cover - optional
        print("[plot] plotly not available; skipping plots.", file=sys.stderr)
        return
    out = Path(outdir or (Path.cwd() / "lacs_output"))
    out.mkdir(parents=True, exist_ok=True)
    for atom, fr in fits.items():
        _plot_atom(atom, fr, out, data_id, method=method, cutoff_k=cutoff_k)


# =============================================================================
# Method-specific fitters (each returns FitResult for one nucleus)
# =============================================================================

def fit_side_rlm_tukey(x: Sequence[float], y: Sequence[float]) -> Tuple[float, float, List[float], List[float]]:
    """Robust linear fit with Tukey biweight via statsmodels RLM.

    Parameters
    ----------
    x, y : sequence[float]
        Data for one side (x ≥ 0 or x < 0).

    Returns
    -------
    (slope, intercept, fitted, residuals) : tuple
        Estimated line parameters and diagnostics.

    Notes
    -----
    Uses IRLS with a redescending ψ; scale estimation uses MAD when available
    (``scale_est='mad'``), with a safe fallback to the default.
    """
    try:
        import statsmodels.api as sm
        from statsmodels.robust.norms import TukeyBiweight
    except Exception as e:
        raise SystemExit("This method requires 'statsmodels' (pip install statsmodels).") from e

    X = sm.add_constant(np.asarray(x, dtype=float))
    yv = np.asarray(y, dtype=float)
    model = sm.RLM(yv, X, M=TukeyBiweight(c=4.685))#optional c=4.685 can be given as tuning paramter for TrkeyBiweight that balances robustness vs. efficiency so the estimator behaves almost like OLS
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
    """Theil–Sen median-of-slopes regression via scikit-learn.

    Parameters
    ----------
    x, y : sequence[float]

    Returns
    -------
    (slope, intercept, fitted, residuals) : tuple

    Notes
    -----
    More robust than OLS with breakdown ≈ 29%, distribution-free; still
    sensitive to extreme x-leverage.
    """
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
    """RANSAC regression with a MAD-derived residual threshold.

    Parameters
    ----------
    x, y : sequence[float]

    Returns
    -------
    (slope, intercept, fitted, residuals) : tuple

    Notes
    -----
    Uses a robust, data-driven residual threshold: baseline OLS residuals -> MAD
    -> σ ≈ 1.4826·MAD -> threshold ≈ 2.5·σ. Compatible with modern sklearn
    (``estimator=...``, numeric ``residual_threshold``).
    """
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
    """Median (τ=0.5) quantile regression via statsmodels.

    Parameters
    ----------
    x, y : sequence[float]

    Returns
    -------
    (slope, intercept, fitted, residuals) : tuple
    """
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
    """Bayesian linear regression with Student‑t noise (PyMC).

    Parameters
    ----------
    x, y : sequence[float]

    Returns
    -------
    (slope, intercept, fitted, residuals, nu, sigma, alpha_draws) : tuple

    Notes
    -----
    Prior defaults (weakly informative on ppm scale):

    - alpha ~ Normal(0, 10)
    - beta  ~ Normal(0, 10)
    - sigma ~ HalfNormal(2)
    - nu    ~ Exponential(1/30) + 1  (df > 1)

    Sampling uses NUTS with reasonable defaults.

    The returned ``alpha_draws`` contains posterior samples of the intercept
    (flattened across chains and draws), which are later used to form offset
    uncertainties.
    """
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

    # Robustly extract arrays from xarray/numpy across PyMC/ArviZ versions
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

    # Flatten posterior draws for alpha (intercept)
    alpha_draws = alpha_arr.reshape(-1)

    fitted = (alpha_mean + beta_mean * x_arr).astype(float).tolist()
    resid = (y_arr - (alpha_mean + beta_mean * x_arr)).astype(float).tolist()
    return beta_mean, alpha_mean, fitted, resid, nu_mean, sigma_mean, alpha_draws


# ---------------------- nucleus-level orchestration ---------------------------

def _fit_atom_by_method(method: str, xvals: List[float], yvals: List[float], tags: List[ResidueKey]) -> FitResult:
    """Fit a single nucleus by splitting x by sign and applying ``method``.

    Parameters
    ----------
    method : {'tukey','theilsen','ransac','quantile','bayes'}
        Robust regression method.
    xvals, yvals : list[float]
        Independent (x) and dependent (y) values per residue.
    tags : list[ResidueKey]
        Residue identifiers in the same order as data.

    Returns
    -------
    FitResult
        Per-side parameters, residuals, and bookkeeping.
    """
    x = np.asarray(xvals, dtype=float); y = np.asarray(yvals, dtype=float)
    pos = x >= 0; neg = ~pos

    def do_side(mask):
        if mask.sum() < 2:
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


def _fit_atom_bayes(xvals: List[float], yvals: List[float], tags: List[ResidueKey]) -> Tuple[FitResult, Dict[str, np.ndarray]]:
    """Bayesian version of :func:`_fit_atom_by_method` that also returns intercept draws.

    Parameters
    ----------
    xvals, yvals : list[float]
        Independent (x) and dependent (y) values per residue.
    tags : list[ResidueKey]
        Residue identifiers in the same order as data.

    Returns
    -------
    (FitResult, dict)
        ``FitResult`` with posterior means, and a dict containing posterior draws
        of the intercept on each side: ``{'pos': alpha_draws_pos, 'neg': alpha_draws_neg}``.
    """
    x = np.asarray(xvals, dtype=float); y = np.asarray(yvals, dtype=float)
    pos = x >= 0; neg = ~pos
    draws: Dict[str, np.ndarray] = {"pos": np.array([], float), "neg": np.array([], float)}

    def do_side(mask, side_key):
        if mask.sum() < 2:
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
             outdir: Optional[Path] = None, plots: bool = True, cutoff_k: float = 2.5) -> Dict[str, Dict]:
    """Run the selected robust method over an NMR-STAR file.

    Parameters
    ----------
    str_file : str
        Path to NMR-STAR file.
    method : {'tukey','theilsen','ransac','quantile','bayes'}
        Robust regression method to apply.
    data_id : str
        Identifier for figure titles/file names.
    rc_model : sequence[str] | str | None
        Random-coil model alias(es) for :class:`RandomCoil`. ``None`` uses the average.
    outdir : Path | None
        Where to save plots; ``None`` defaults to ``./lacs_output``.
    plots : bool
        Whether to generate Plotly plots.
    cutoff_k : float
        Outlier cutoff multiplier (k in |r|/(k·MAD)).

    Returns
    -------
    dict
        Structured report keyed by chemical-shift list id.
    """
    cs = read_star(str_file)
    results: Dict[str, Dict] = {}
    for list_id, resmap in cs.items():
        d, tags = compute_deltas(resmap, rc_model)
        fits: Dict[str, FitResult] = {}
        alpha_samples: Dict[str, Dict[str, np.ndarray]] = {}

        for atom, xkey in [('ca','ca'), ('cb','ca'), ('c','x_for_c'), ('n','x_for_n'), ('ha','x_for_ha')]:
            yvals = d[atom]
            if atom in {'ca','cb'}:
                xvals = [a - b for a, b in zip(d['ca'], d['cb'])]
                tg = tags[atom]
            else:
                xvals = d[xkey]; tg = tags[atom]
            if len(yvals) >= 2 and len(xvals) == len(yvals):
                if method == "bayes":
                    fr, draws = _fit_atom_bayes(xvals, yvals, tg)
                    alpha_samples[atom] = draws
                else:
                    fr = _fit_atom_by_method(method, xvals, yvals, tg)
                fits[atom] = fr

        if fits:
            maybe_plot_all(fits, outdir, data_id, method, plots, cutoff_k=cutoff_k)
            if method == "bayes":
                results[list_id] = collect_and_report_bayes(fits, alpha_samples, cutoff_k=cutoff_k)
            else:
                results[list_id] = collect_and_report(fits, cutoff_k=cutoff_k)
    return results


# =============================================================================
# CLI
# =============================================================================

def _default_json_path(outdir: Optional[Path], data_id: str, method: str) -> Path:
    """Compute a default JSON output filepath.

    Parameters
    ----------
    outdir : Path | None
        Base directory (defaults to CWD if ``None``).
    data_id : str
        Identifier for the dataset/entry.
    method : str
        Method name used in the run.

    Returns
    -------
    Path
        Full path to ``<data_id>_<method>.json``.
    """
    base = Path.cwd() if outdir is None else Path(outdir)
    return base / f"{data_id}_{method}.json"


def main(argv=None) -> None:
    """Command-line entry point.

    Parameters
    ----------
    argv : list[str] | None
        Argument vector; if ``None``, uses ``sys.argv[1:]``.

    Side Effects
    ------------
    - Generates optional Plotly plots (HTML and PDF if ``kaleido`` is installed).
    - Writes a JSON report to ``--json-out`` or to a sensible default.

    Examples
    --------
    .. code-block:: bash

        python lacs_unified.py entry.str --method theilsen --data-id demo --out figs --json-out demo_theilsen.json
    """
    import argparse

    p = argparse.ArgumentParser(description="LACS unified CLI for robust linear fits.")
    p.add_argument("star_file", help="Path to NMR-STAR .str file")
    p.add_argument("--method", required=True,
                   choices=["tukey","theilsen","ransac","quantile","bayes"],
                   help="Robust regression method to use")
    p.add_argument("--data-id", default="LACS")
    p.add_argument("--rc-model", nargs="*", default=None,
                   help="Random-coil model alias(es), e.g. wis wan; omit for average of all")
    p.add_argument("--out", type=Path, default=None, help="Output directory for plots")
    p.add_argument("--no-plots", action="store_true", help="Disable plotting entirely")
    p.add_argument("--cutoff-k", type=float, default=5.0, help="Outlier cutoff multiplier (default 5.0)")
    p.add_argument("--json-out", type=Path, default=None, help="Where to write JSON (defaults to <data_id>_<method>.json)")

    args = p.parse_args(argv)
    rc_model = None if args.rc_model == [] else (args.rc_model if args.rc_model is not None else None)
    report = run_lacs(args.star_file, method=args.method, data_id=args.data_id,
                      rc_model=rc_model, outdir=args.out, plots=not args.no_plots, cutoff_k=args.cutoff_k)

    json_path = args.json_out or _default_json_path(args.out, args.data_id, args.method)
    json_path.parent.mkdir(parents=True, exist_ok=True)
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2)
    print(str(json_path))


if __name__ == "__main__":
    main()
