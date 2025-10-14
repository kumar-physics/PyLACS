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

from typing import Dict, List, Optional, Sequence, Tuple
from typing import Dict, Any, List

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

try:
    from pylacs.apply_lacs_correction import apply_selected_offsets_and_note
except Exception as e:
    raise SystemExit("Couldn'timport apply_lacs_correction.py. Ensure it is on PYTHONPATH or in the same folder.") from e
from datetime import datetime, timezone
try:
    # Python 3.11+
    from datetime import UTC
except ImportError:
    # Older versions
    UTC = timezone.utc
from pathlib import Path
from typing import Dict, Any, Optional, Sequence
import json
import numpy as np


def _star_tag_value(value: Any):
    """Serializer for SAVEFRAME TAGS: None for empty -> '.' in STAR."""
    try:
        from pathlib import Path
    except Exception:
        Path = None  # type: ignore
    import numpy as np
    if value is None:
        return None
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, (list, tuple, set)):
        if len(value) == 0:
            return None
        value = " ".join(map(str, value))
    if Path and isinstance(value, Path):
        value = str(value)
    s = str(value).strip()
    return None if s == "" else s

def _star_loop_cell(value: Any) -> str:
    """Serializer for LOOP CELLS: '.' for missing/empty."""
    try:
        from pathlib import Path
    except Exception:
        Path = None  # type: ignore
    import numpy as np
    if value is None:
        return "."
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, (list, tuple, set)):
        if len(value) == 0:
            return "."
        value = " ".join(map(str, value))
    if Path and isinstance(value, Path):
        value = str(value)
    s = str(value).strip()
    return "." if s == "" else s

def _make_loop(p, category: str, tag_names: List[str]):
    """Create a pynmrstar Loop robustly across versions."""
    for kw in ("tag_names", "tags"):
        try:
            return p.Loop.from_scratch(category=category, **{kw: tag_names})
        except TypeError:
            pass
    loop = p.Loop.from_scratch(category=category)
    if hasattr(loop, "set_tag_names"):
        loop.set_tag_names(tag_names)
    elif hasattr(loop, "add_tag_name"):
        for nm in tag_names:
            loop.add_tag_name(nm)
    elif hasattr(loop, "add_tag"):
        for nm in tag_names:
            loop.add_tag(nm)
    return loop

def _flatten_meta(d: Dict[str, Any], parent: str = "") -> List[tuple[str, str]]:
    """Flatten a small dict (e.g., meta) into (key_path, value_str)."""
    rows: List[tuple[str, str]] = []
    for k, v in (d or {}).items():
        key = f"{parent}.{k}" if parent else k
        if isinstance(v, dict):
            rows.extend(_flatten_meta(v, key))
        else:
            rows.append((key, _star_loop_cell(v)))
    return rows

def write_star_pynmrstar(report: Dict[str, Any], filepath: str, params: Dict[str, Any]) -> None:
    """
      :param report: output dictonary frm run_lacs
      :param filepath: outut STAR file path
      :param params: input parameters of run_lacs
      :return: None (writes outout file)

      Write STAR with requested prefixes:

      - Metadata:   _LACS_metadata.
      - Offsets:    _LACS_offsets.     (includes pos/neg and optional Bayes overall row)
      - Fit data:   _LACS_fitdata.    (x, y, residual, flag, prob)
      - Per-list meta (from JSON 'meta') also under _LACS_metadata. loop

    """
    import pynmrstar as p

    entry = p.Entry.from_scratch("pylacs_lacs")

    # --------------------- global metadata (run params) -----------------------
    sf_meta = p.Saveframe.from_scratch("save_LACS_metadata")
    sf_meta.add_tag("_LACS_metadata.Generated_by", "PyLACS")
    sf_meta.add_tag("_LACS_metadata.Generation_date", datetime.now(UTC).isoformat())
    for k, v in (params or {}).items():
        sf_meta.add_tag(f"_LACS_metadata.{k}", _star_tag_value(v))
    entry.add_saveframe(sf_meta)

    # ---------------------- per-list metadata (from JSON) ---------------------
    # Mirrors JSON 'meta' dict (e.g., cutoff_k, used_side, bayes=True, etc.)
    meta_tags = [
        "_LACS_metadata.List_ID",
        "_LACS_metadata.Key",
        "_LACS_metadata.Value",
    ]
    sf_meta_pl = p.Saveframe.from_scratch("save_LACS_metadata_perlist")
    sf_meta_pl.category   = "LACS_metadata"   # loops only → set category & prefix
    sf_meta_pl.tag_prefix = "_LACS_metadata"
    loop_meta_pl = _make_loop(p, category="LACS_metadata", tag_names=meta_tags)

    meta_rows: List[List[str]] = []
    for list_id, rep in (report or {}).items():
        flat = _flatten_meta(rep.get("meta", {}))
        for key_path, val in flat:
            meta_rows.append([_star_loop_cell(list_id), _star_loop_cell(key_path), val])
    if meta_rows:
        loop_meta_pl.add_data(meta_rows)
    sf_meta_pl.add_loop(loop_meta_pl)
    entry.add_saveframe(sf_meta_pl)
    # ---------------------- Offsets (pos/neg + Bayes overall) ----------------
    off_tags = [
        "_LACS_offsets.List_ID",
        "_LACS_offsets.Atom",
        "_LACS_offsets.Side",  # 'pos' | 'neg' | 'overall'
        "_LACS_offsets.Value",  # pos/neg split value or Bayes mean for 'overall'
        "_LACS_offsets.Bayes_mean",
        "_LACS_offsets.Bayes_ci95_lo",
        "_LACS_offsets.Bayes_ci95_hi",
        "_LACS_offsets.Bayes_sd",
    ]
    sf_off = p.Saveframe.from_scratch("save_LACS_offsets")
    sf_off.category = "LACS_offsets"
    sf_off.tag_prefix = "_LACS_offsets"
    loop_off = _make_loop(p, category="LACS_offsets", tag_names=off_tags)

    off_rows: List[List[str]] = []
    for list_id, rep in (report or {}).items():
        split = rep.get("offsets_split", {})  # deterministic per-side offsets
        bayes_overall = rep.get("offsets_bayes", {})  # overall Bayes stats
        # Accept either name for per-side Bayes stats:
        bayes_sides = (rep.get("offsets_bayes_sides") or
                       rep.get("offsets_bayes_split") or {})  # per-side Bayes stats

        # 1) pos/neg rows (fill Bayes_* from per-side stats if present)
        for atom, sides in (split or {}).items():
            for side in ("pos", "neg"):
                val = None if sides is None else sides.get(side)
                if val is None:
                    continue
                b = ((bayes_sides.get(atom) or {}).get(side) or {})
                ci = b.get("ci95") or [None, None]
                off_rows.append([
                    _star_loop_cell(list_id),
                    _star_loop_cell(str(atom).upper()),
                    _star_loop_cell(side),
                    _star_loop_cell(val),
                    _star_loop_cell(b.get("mean")),
                    _star_loop_cell(ci[0]),
                    _star_loop_cell(ci[1]),
                    _star_loop_cell(b.get("sd")),
                ])

        # 2) Bayes "overall" row
        for atom, stats in (bayes_overall or {}).items():
            mean = stats.get("mean")
            ci95 = stats.get("ci95") or [None, None]
            sd = stats.get("sd")
            off_rows.append([
                _star_loop_cell(list_id),
                _star_loop_cell(str(atom).upper()),
                _star_loop_cell("overall"),
                _star_loop_cell(mean if mean is not None else "."),  # duplicate mean in Value
                _star_loop_cell(mean),
                _star_loop_cell(ci95[0]),
                _star_loop_cell(ci95[1]),
                _star_loop_cell(sd),
            ])

    if off_rows:
        loop_off.add_data(off_rows)
    sf_off.add_loop(loop_off)
    entry.add_saveframe(sf_off)

    # -------------------- Fit data (x,y,residual,flag,prob) -------------------
    fit_tags = [
        "_LACS_fitdata.List_ID",
        "_LACS_fitdata.Atom",
        "_LACS_fitdata.Entity_ID",
        "_LACS_fitdata.Entity_assembly_ID",
        "_LACS_fitdata.Comp_index_ID",
        "_LACS_fitdata.Comp_ID",
        "_LACS_fitdata.X",
        "_LACS_fitdata.Y",
        "_LACS_fitdata.Residual",
        "_LACS_fitdata.Flag",
        "_LACS_fitdata.Prob",
    ]
    sf_fit = p.Saveframe.from_scratch("save_LACS_fitdata")
    sf_fit.category   = "LACS_fitdata"
    sf_fit.tag_prefix = "_LACS_fitdata"
    loop_fit = _make_loop(p, category="LACS_fitdata", tag_names=fit_tags)

    fit_rows: List[List[str]] = []
    for list_id, rep in (report or {}).items():
        outliers = rep.get("outliers", {})
        for atom, rows in (outliers or {}).items():
            for row in rows:
                rk = row.get("residue_key", {})
                fit_rows.append([
                    _star_loop_cell(list_id),
                    _star_loop_cell(str(atom).upper()),
                    _star_loop_cell(rk.get("entity")),
                    _star_loop_cell(rk.get("assembly")),
                    _star_loop_cell(rk.get("index")),
                    _star_loop_cell(rk.get("comp")),
                    _star_loop_cell(row.get("x")),
                    _star_loop_cell(row.get("y")),
                    _star_loop_cell(row.get("residual")),
                    _star_loop_cell(row.get("flag")),
                    _star_loop_cell(row.get("prob")),
                ])
    if fit_rows:
        loop_fit.add_data(fit_rows)
    sf_fit.add_loop(loop_fit)
    entry.add_saveframe(sf_fit)

    # ----------------------- write -----------------------
    entry.write_to_file(filepath, skip_empty_tags=True, skip_empty_loops=True)

def write_report(
    report: Dict[str, Any],
    base_path: Path,
    write_format: str = "both",
    params_for_star: Optional[Dict[str, Any]] = None,
    json_out: Optional[Path] = None,
    star_out: Optional[Path] = None,
) -> Dict[str, Path]:
    """
    Write report to disk as JSON, STAR, or both.
    Returns a dict of the actual paths written, e.g. {'json': Path(...), 'star': Path(...)}.

    :param report: Output dictionary from run_lacs
    :param base_path: Output path
    :param write_format: 'json' or 'star' or 'both' default 'both'
    :param params_for_star: Input parameters of run_lacs
    :param json_out: JSON output file path
    :param star_out: STAR output file path
    :return: Returns a dict of the actual paths written, e.g. {'json': Path(...), 'star': Path(...)}.
    """
    write_format = (write_format or "json").lower()
    out_paths: Dict[str, Path] = {}

    if write_format in ("json", "both"):
        jp = Path(json_out) if json_out else base_path.with_suffix(".json")
        jp.parent.mkdir(parents=True, exist_ok=True)
        with open(jp, "w", encoding="utf-8") as f:
            json.dump(report, f, indent=2)
        out_paths["json"] = jp

    if write_format in ("star", "both"):
        sp = Path(star_out) if star_out else base_path.with_suffix(".star")
        write_star_pynmrstar(report, sp, params_for_star or {})
        out_paths["star"] = sp

    if not out_paths:
        raise ValueError("write_format must be one of {'json','star','both'}.")

    return out_paths

def apply_corrections_from_report(
    report: Dict[str, Any],
    input_star: Path,
    output_star: Path,
    data_id: str,
    atoms: Sequence[str] = ("CA", "CB", "C", "N"),
    release_author: str = "BMRB",
) -> Dict[str, int]:
    """
    Apply LACS offsets to an input NMR-STAR file for each list_id in report, and write a corrected file.

    :param report: LACS report returned by run_lacs(...).
    :param input_star: Source .str/.star file.
    :param output_star: Destination corrected .str/.star file.
    :param data_id: Identifier used in the Release.Detail note.
    :param atoms: Subset of atoms to correct (default CA, CB, C, N).
    :param release_author: Author recorded in Release loop.
    :return: Counts across the last applied list_id; includes per-atom keys and 'total'.
        (If multiple list_ids exist, counts reflect the most recent application.)
    """
    if apply_selected_offsets_and_note is None:
        raise RuntimeError("apply_selected_offsets_and_note not available. Ensure apply_lacs_correction.py is importable.")

    atoms_use = _normalize_atoms(atoms)

    # Ensure string paths for pynmrstar compatibility
    input_star_s = str(input_star)
    output_star_s = str(output_star)

    Path(output_star_s).parent.mkdir(parents=True, exist_ok=True)
    if not Path(input_star_s).is_file():
        raise FileNotFoundError(f"Input STAR file not found: {input_star_s}")
    if Path(output_star_s).suffix.lower() not in {".str", ".star"}:
        output_star_s = str(Path(output_star_s).with_suffix(".str"))

    counts_last: Dict[str, int] = {}
    # Iterate in the order present in the report
    for list_id in report:
        offsets_uc = _extract_offsets_for_list(report, list_id)
        parts = [f"{a}={offsets_uc[a]:+g}" for a in atoms_use]
        details = (
            f"LACS correction applied to list_id {list_id} ({data_id}): "
            + (", ".join(parts) if parts else "none")
            + ". Source: LACS."
        )
        counts_last = apply_selected_offsets_and_note(
            input_path=input_star_s,
            output_path=output_star_s,
            list_id=int(list_id),
            offsets=offsets_uc,
            atoms=atoms_use,
            release_author=release_author,
            release_details=details,
        )
    return counts_last


def _extract_offsets_for_list(report: Dict, list_id: int) -> Dict[str, float]:
    """
    Given a full LACS report (dict loaded from JSON or returned by run_lacs),
    extract offsets for the specified list_id.
    Returns keys in UPPER CASE: {'CA','CB','C','N'} (missing → 0.0).


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


    # Map lower-case keys used by LACS ('ca','cb','c','n') → upper-case
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
    """
    this
    :param atoms:
    :return:
    """

    allowed = {"CA", "CB", "C", "N"}
    if not atoms:
        return ["CA", "CB", "C", "N"]
    atoms_up = [a.upper() for a in atoms]
    return [a for a in atoms_up if a in allowed]

# ---------------------------------------------------------------------------
# NEW PUBLIC WRAPPER: apply correction from an existing JSON report
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

    :param str_file: Path to input NMR-STAR file.
    :param data_id: Data ID to place in the Release row.
    :param lacs_output: Path to LACS JSON output file.
    :param output_corrected: Where to write corrected file.
    :param list_id: Chemical-shift list ID to modify.
    :param atoms: Subset of atoms to which offsets should be applied (defaults to all four).
    :param release_author: Author to place in the Release row (default 'BMRB').
    :return: Counts dictionary with per-atom and total updates.
    """
    if apply_selected_offsets_and_note is None:
        raise RuntimeError("apply_selected_offsets_and_note not available. Ensure apply_lacs_correction.py is on PYTHONPATH.")

    with open(lacs_output, "r", encoding="utf-8") as f:
        report = json.load(f)

    offsets_uc = _extract_offsets_for_list(report, list_id)
    atoms_use = _normalize_atoms(atoms)

    # Compose a concise details string for the Release.Detail cell
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
#Place holder for scale factors, if available in the future
CA_SCALE = {'ILE': 4.8, 'GLN': 4.5, 'GLY': 1.7, 'GLU': 4.5, 'CYS': 6.6, 'ASP': 4.2, 'SER': 4.1, 'LYS': 4.6, 'PRO': 3.1, 'ASN': 3.6, 'VAL': 5.8, 'THR': 4.9, 'HIS': 4.0, 'TRP': 4.0, 'PHE': 5.3, 'ALA': 5.1, 'MET': 4.2, 'LEU': 4.3, 'ARG': 4.6, 'TYR': 4.7}
CB_SCALE = {'ILE': 2.4, 'GLN': 3.5, 'GLY': 1.0, 'GLU': 3.4, 'CYS': 3.6, 'ASP': 1.9, 'SER': 2.55, 'LYS': 2.4, 'PRO': 1.0, 'ASN': 1.7, 'VAL': 2.7, 'THR': 1.5, 'HIS': 3.4, 'TRP': 1.9, 'PHE': 2.5, 'ALA': 4.2, 'MET': 2.8, 'LEU': 2.65, 'ARG': 2.7, 'TYR': 3.1}
POS_SCALE = {'ILE': 4.83, 'GLN': 3.82, 'GLY': 1.0, 'GLU': 3.24, 'CYS': 4.35, 'ASP': 3.42, 'SER': 3.98, 'LYS': 3.42, 'PRO': 3.12, 'ASN': 2.94, 'VAL': 5.19, 'THR': 5.56, 'HIS': 4.26, 'TRP': 3.09, 'PHE': 4.04, 'ALA': 4.01, 'MET': 3.44, 'LEU': 3.47, 'ARG': 3.65, 'TYR': 3.84}
NEG_SCALE = {'ILE': 2.37, 'GLN': 4.18, 'GLY': 1.0, 'GLU': 4.66, 'CYS': 5.85, 'ASP': 2.68, 'SER': 2.67, 'LYS': 3.58, 'PRO': 0.98, 'ASN': 2.36, 'VAL': 3.31, 'THR': 0.84, 'HIS': 3.14, 'TRP': 2.81, 'PHE': 3.76, 'ALA': 5.29, 'MET': 3.56, 'LEU': 3.48, 'ARG': 3.65, 'TYR': 3.96}

def _tag_to_label(tag: tuple) -> str:
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

    :param slope_pos: Positive slope.
    :param intercept_pos: Positive intercept.
    :param fitted_pos: Fitted values for positive side.
    :param resid_pos: Residuals for positive side.
    :param x_pos: X values for positive side.
    :param y_pos: Y values for positive side.
    :param tags_pos: Residue tags for positive side.
    :param slope_neg: Negative slope.
    :param intercept_neg: Negative intercept.
    :param fitted_neg: Fitted values for negative side.
    :param resid_neg: Residuals for negative side.
    :param x_neg: X values for negative side.
    :param y_neg: Y values for negative side.
    :param tags_neg: Residue tags for negative side.


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

    :param file_name: Path to NMR-STAR file.
    :return: Chemical-shift map for each list id as a dict of dicts.
        Keys are residue tags, values are dicts of atom names to chemical shifts.

    Notes:
    - Only rows in the ``Atom_chem_shift`` category are parsed. Non-numeric
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

    :param resmap: Chemical-shift map for each list id as a dict of dicts.
        Keys are residue tags, values are dicts of atom names to chemical shifts.
    :param rc_model: Random-coil model to use.
    :return: Tuple of dicts of lists and dicts of lists of residue tags.
        Keys are atom names, values are lists of Δδ values and residue tags.
    :raises ValueError: If the random-coil model is not recognized.
    :raises KeyError: If the chemical-shift map is empty.

    Notes:

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

def compute_deltas2(resmap: Dict[ResidueKey, Dict[str, float]], rc_model: Optional[Sequence[str] | str]) -> Tuple[Dict[str, List[float]], Dict[str, List[ResidueKey]]]:
    """Compute Δδ against random-coil and assemble arrays for each nucleus.

    :param resmap: Chemical-shift map for each list id as a dict of dicts.
        Keys are residue tags, values are dicts of atom names to chemical shifts.
    :param rc_model: Random-coil model to use.
    :return: Tuple of dicts of lists and dicts of lists of residue tags.
        Keys are atom names, values are lists of Δδ values and residue tags.
    :raises ValueError: If the random-coil model is not recognized.
    :raises KeyError: If the chemical-shift map is empty.

    Notes:

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

            if ca is not None and cb is not None and abs(ca)<15  and 15 > cb > -25:
                x = ca - cb
                # if x1 >=0 :
                #     x=x1/POS_SCALE[comp_id]
                # else:
                #     x=x1/NEG_SCALE[comp_id]


                d['ca'].append(ca); tags['ca'].append(residue)
                d['cb'].append(cb); tags['cb'].append(residue)
                if c is not None and abs(c)<12.0:
                    d['c'].append(c); d['x_for_c'].append(x); tags['c'].append(residue)
                if n is not None and abs(n)<20.0:
                    d['n'].append(n); d['x_for_n'].append(x); tags['n'].append(residue)
                if ha is not None and abs(ha)<4.0:
                    d['ha'].append(ha); d['x_for_ha'].append(x); tags['ha'].append(residue)
        except Exception:
            continue
    return d, tags

def mad(arr: np.ndarray) -> float:
    """Median absolute deviation (MAD).

    :param arr: Array of values.
    :return: Median absolute deviation.
    :raises ValueError: If the array is empty.


    """
    med = np.median(arr)
    m = np.median(np.abs(arr - med))
    return float(m) if m > 0 else 1.0


def logistic_prob(z: np.ndarray, slope: float = 6.0, hinge: float = 1.0) -> np.ndarray:
    """Logistic map from a nonnegative score to a soft probability in [0, 1].

    :param z: Array of scores.
    :param slope: Slope of the logistic map (default 6.0).
    :param hinge: Hinge of the logistic map (default 1.0).
    :return: Array of soft probabilities.
    :raises ValueError: If the array is empty.
    """
    return 1.0 / (1.0 + np.exp(-slope * (z - hinge)))


def outlier_stats(residuals: np.ndarray, scale: Optional[float] = None, cutoff_k: float = 5.0) -> Tuple[List[int], List[float]]:
    """Compute 0/1 outlier flags and smooth probabilities from residuals.

    :param residuals: Array of residuals.
    :param scale: Scale factor for outlier detection (default 5.0).
    :param cutoff_k: Cutoff multiplier for outlier detection (default 5.0).
    :return: Tuple of lists of 0/1 flags and soft probabilities.
    :raises ValueError: If the array is empty.
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

    If both sides were fit (non-empty), offset = -average(intercepts).
    If only one side was fit, offset = -that side's intercept.

    :param fits: Dictionary of fit results.
    :param cutoff_k: Outlier cutoff multiplier used for residual-based probabilities.
    :return: Dictionary of offsets and outlier lists.
    """
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
        # >>> NEW: x/y aligned to tags across both sides
        x_all = (fr.x_pos + fr.x_neg)
        y_all = (fr.y_pos + fr.y_neg)
        # <<< NEW
        flags, probs = outlier_stats(resid_all, cutoff_k=cutoff_k)
        out: List[Dict[str, object]] = []
        for tag, f, p, r, xv, yv in zip(tags_all, flags, probs, resid_all.tolist(), x_all, y_all):
            out.append({
                "residue_key": {"entity": tag[0], "assembly": tag[1], "index": tag[2], "comp": tag[3]},
                "residual": round(float(r), 4),
                "x": round(float(xv), 4),  # <<< NEW
                "y": round(float(yv), 4),  # <<< NEW
                "flag": int(f),
                "prob": round(float(p), 4),
            })
        report["outliers"][atom] = out
    return report


def collect_and_report_bayes(
    fits: Dict[str, FitResult],
    alpha_samples: Dict[str, Dict[str, np.ndarray]],
    cutoff_k: float = 5.0
) -> Dict[str, Dict]:
    """
    Assemble a full Bayes report.

    :param fits: Dictionary of fit results.
    :param cutoff_k: Outlier cutoff multiplier used for residual-based probabilities.
    :return: Dictionary of offsets and outlier lists.

    Returns a dict with keys:
      - offsets_bayes:        {atom: {"mean": float, "ci95": [lo, hi], "sd": float|None}}
      - offsets_bayes_sides:  {atom: {"pos": {...}, "neg": {...}}}
      - offsets_split:        {atom: {"pos": float|None, "neg": float|None}}
      - outliers:             {atom: [ {residue_key: {...}, x, y, residual, flag, prob}, ... ]}
      - meta:                 {"cutoff_k": float, "bayes": True, "used_side": {atom: "both"|"pos"|"neg"|"none"}}
    """
    import numpy as np

    offsets_bayes: Dict[str, Dict[str, Any]] = {}
    offsets_bayes_sides: Dict[str, Dict[str, Dict[str, Optional[float]]]] = {}
    offsets_split: Dict[str, Dict[str, Optional[float]]] = {}
    used_side: Dict[str, str] = {}
    outliers_map: Dict[str, List[Dict[str, object]]] = {}

    for atom, fr in fits.items():
        # Which sides have data?
        has_pos = len(fr.x_pos) > 0
        has_neg = len(fr.x_neg) > 0
        used_side[atom] = "both" if (has_pos and has_neg) else ("pos" if has_pos else ("neg" if has_neg else "none"))

        # Deterministic per-side offsets from fitted intercepts
        offsets_split[atom] = {
            "pos": round(float(-fr.intercept_pos), 3) if has_pos else None,
            "neg": round(float(-fr.intercept_neg), 3) if has_neg else None,
        }

        # Posterior draws for offset (= -alpha) per side
        pos_draws = -alpha_samples.get(atom, {}).get("pos", np.array([], float))
        neg_draws = -alpha_samples.get(atom, {}).get("neg", np.array([], float))

        # Overall draws: combine sides when both present
        if pos_draws.size > 0 and neg_draws.size > 0:
            n = min(pos_draws.size, neg_draws.size)
            overall_draws = 0.5 * (pos_draws[:n] + neg_draws[:n])
        elif pos_draws.size > 0:
            overall_draws = pos_draws
        elif neg_draws.size > 0:
            overall_draws = neg_draws
        else:
            overall_draws = np.array([], float)

        # Overall Bayes stats (fallback to intercepts if no draws)
        if overall_draws.size:
            mean = float(np.mean(overall_draws))
            lo, hi, sd = _ci_and_sd(overall_draws, level=0.95)  # returns (lo, hi, sd)
        else:
            mean = float(-0.5 * (fr.intercept_pos + fr.intercept_neg))
            lo = hi = sd = float("nan")

        offsets_bayes[atom] = {
            "mean": round(mean, 4),
            "ci95": [None if np.isnan(lo) else round(lo, 4),
                     None if np.isnan(hi) else round(hi, 4)],
            "sd":   None if np.isnan(sd) else round(sd, 4),
        }

        # Per-side Bayes stats
        def _side_stats(arr: np.ndarray) -> Dict[str, Optional[float] | List[Optional[float]]]:
            if arr.size == 0:
                return {"mean": None, "ci95": [None, None], "sd": None}
            m = float(np.mean(arr))
            lo_s, hi_s, sd_s = _ci_and_sd(arr, level=0.95)
            return {
                "mean": round(m, 4),
                "ci95": [None if np.isnan(lo_s) else round(lo_s, 4),
                         None if np.isnan(hi_s) else round(hi_s, 4)],
                "sd": None if np.isnan(sd_s) else round(sd_s, 4),
            }

        offsets_bayes_sides[atom] = {
            "pos": _side_stats(pos_draws),
            "neg": _side_stats(neg_draws),
        }

        # Outliers (with x/y attached)
        resid_all = np.array(fr.resid_pos + fr.resid_neg, dtype=float)
        tags_all  = fr.tags_pos + fr.tags_neg
        x_all     = fr.x_pos + fr.x_neg
        y_all     = fr.y_pos + fr.y_neg

        flags, probs = outlier_stats(resid_all, cutoff_k=cutoff_k)
        rows: List[Dict[str, object]] = []
        for tag, f, p, r, xv, yv in zip(tags_all, flags, probs, resid_all.tolist(), x_all, y_all):
            rows.append({
                "residue_key": {"entity": tag[0], "assembly": tag[1], "index": tag[2], "comp": tag[3]},
                "x": round(float(xv), 6),
                "y": round(float(yv), 6),
                "residual": round(float(r), 4),
                "flag": int(f),
                "prob": round(float(p), 4),
            })
        outliers_map[atom] = rows

    return {
        "offsets_bayes": offsets_bayes,
        "offsets_bayes_sides": offsets_bayes_sides,
        "offsets_split": offsets_split,
        "outliers": outliers_map,
        "meta": {"cutoff_k": cutoff_k, "bayes": True, "used_side": used_side},
    }

# =============================================================================
# Plotting (optional): scatter with outlier marks and probability bars
# =============================================================================


def _plot_atom(atom: str, fr: FitResult, outdir: Path, data_id: str, list_id, method: str, cutoff_k: float = 5.0) -> None:
    """Render interactive plots for a single nucleus with residue labels.

   :param atom: Atom name.
   :param fr: Fit result.
   :param outdir: Output directory.
   :param data_id: Data identifier.
   :param list_id: List identifier.
   :param method: Method name.
   :param cutoff_k: Outlier cutoff.
    """
    if not PLOTLY_AVAILABLE:
        return

    x_pos, y_pos = np.asarray(fr.x_pos, float), np.asarray(fr.y_pos, float)
    x_neg, y_neg = np.asarray(fr.x_neg, float), np.asarray(fr.y_neg, float)

    # --- NEW: if absolutely no points, skip plotting entirely
    if (x_pos.size + x_neg.size) == 0:
        return
    # Pre-compute labels like "56H"
    labels_pos = [_tag_to_label(t) for t in fr.tags_pos]
    labels_neg = [_tag_to_label(t) for t in fr.tags_neg]

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


def maybe_plot_all(fits: Dict[str, FitResult], outdir: Optional[Path], data_id: str, method: str, enable_plots: bool, list_id,cutoff_k: float = 5.0) -> None:
    """Render all per-nucleus plots if plotting is enabled and supported.

    :param fits: Dictionary of FitResult objects for each nucleus.
    :param outdir: Output directory for plots.  If ``None``, ``./lacs_output`` is used.
    :param data_id: Data identifier.
    :param method: Method identifier.
    :param enable_plots: Whether to enable plotting. If ``False`` or Plotly is unavailable, plotting is skipped.
    :param list_id: List identifier.
    :param cutoff_k: Cutoff for outlier detection.


    """
    if not enable_plots or not fits:
        return
    if not PLOTLY_AVAILABLE:  # pragma: no cover - optional
        print("[plot] plotly not available; skipping plots.", file=sys.stderr)
        return
    out = Path(outdir or (Path.cwd() / "lacs_output"))
    out.mkdir(parents=True, exist_ok=True)
    for atom, fr in fits.items():
        _plot_atom(atom, fr, out, data_id, list_id,method=method, cutoff_k=cutoff_k)


# =============================================================================
# Method-specific fitters (each returns FitResult for one nucleus)
# =============================================================================

def fit_side_rlm_tukey(x: Sequence[float], y: Sequence[float]) -> Tuple[float, float, List[float], List[float]]:
    """Robust linear fit with Tukey biweight via statsmodels RLM.

    :param x: sequence[float]
    :param y: sequence[float]
    :return: (slope, intercept, fitted, residuals) : tuple

    Notes:
    - Uses IRLS with a redescending ψ; scale estimation uses MAD when available
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

    :param x: sequence[float]
    :param y: sequence[float]
    :return: (slope, intercept, fitted, residuals) : tuple

    Notes:
    - More robust than OLS with breakdown ≈ 29%, distribution-free; still
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

    :param x: sequence[float]
    :param y: sequence[float]
    :return: (slope, intercept, fitted, residuals) : tuple

    Notes:
    - Uses a robust, data-driven residual threshold: baseline OLS residuals -> MAD
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

    :param x: sequence[float]
    :param y: sequence[float]
    :return: (slope, intercept, fitted, residuals) : tuple

    Notes:
    - More robust than OLS with breakdown ≈ 50%, distribution-free; still
    sensitive to extreme x-leverage.
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

    :param x: sequence[float]
    :param y: sequence[float]
    :return: (slope, intercept, fitted, residuals) : tuple

    Notes:

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

def _fit_atom_by_method(method: str, xvals: List[float], yvals: List[float], tags: List[ResidueKey],
                        min_per_side: int = MIN_PER_SIDE_DEFAULT) -> FitResult:
    """Fit a single nucleus by splitting x by sign and applying ``method``.

    Only fit a side if it contains at least ``min_per_side`` points.
    """
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
    """Bayesian version of :func:`_fit_atom_by_method` that also returns intercept draws."""
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


def run_lacs2(cs : Dict[str, Dict[ResidueKey, Dict[str, float]]] = {}, method: str='tukey', data_id: str ='BMRB', rc_model: Optional[Sequence[str] | str] = None,
             outdir: Optional[Path] = None, plots: bool = True, cutoff_k: float = 5.0,
             min_per_side: int = 5,
             write_format: str = "json",  # {'json','star','both'}
             json_out: Optional[Path] = None,
             star_out: Optional[Path] = None,
             params_for_star: Optional[Dict[str, Any]] = None,  # if None, we'll build one internally
             apply_corrections: bool = False,
             correction_atoms: Sequence[str] = ("CA", "CB", "C", "N"),
             release_author: str = "BMRB",
             output_corrected: Optional[Path] = None) -> Dict[str, Dict]:
    """Run the selected robust method over an NMR-STAR file.

    :param str_file: Path to NMR-STAR file.
    :param method: Robust regression method to use.
    :param data_id: Identifier for the dataset/entry.
    :param rc_model: Random-coil model alias(es), e.g. wis wan; omit for average of all.
    :param outdir: Output directory for plots.``None`` defaults to ``./lacs_output``.
    :param plots: Whether to generate plots.
    :param cutoff_k: Outlier cutoff multiplier.(k in |r|/(k·MAD)).
    :param min_per_side: Minimum number of points required on each sign side.
    :return: Dictionary of results, keyed by list ID.

    """

    #cs = read_star(str_file)
    results: Dict[str, Dict] = {}
    for list_id, resmap in cs.items():
        d, tags = compute_deltas2(resmap, rc_model)
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
                    fr, draws = _fit_atom_bayes(xvals, yvals, tg, min_per_side=min_per_side)
                    alpha_samples[atom] = draws
                else:
                    fr = _fit_atom_by_method(method, xvals, yvals, tg, min_per_side=min_per_side)
                fits[atom] = fr

        if fits:
            maybe_plot_all(fits, outdir, data_id, method, plots, list_id,cutoff_k=cutoff_k)
            if method == "bayes":
                results[list_id] = collect_and_report_bayes(fits, alpha_samples, cutoff_k=cutoff_k)
            else:
                results[list_id] = collect_and_report(fits, cutoff_k=cutoff_k)
    # ---------- Writing (JSON / STAR / both) ----------
    base_dir = Path.cwd() if outdir is None else Path(outdir)
    base_dir.mkdir(parents=True, exist_ok=True)
    base_path = base_dir / f"{data_id}_{method}"
    report = results
    # If caller didn't pass STAR metadata, build a minimal one from our args
    if params_for_star is None:
        params_for_star = dict(
            str_file=str_file,
            method=method,
            data_id=data_id,
            rc_model=rc_model if rc_model is not None else "",
            outdir=str(base_dir),
            plots=plots,
            cutoff_k=cutoff_k,
            min_per_side=min_per_side,
        )

    written = write_report(
        report=report,
        base_path=base_path,
        write_format=write_format,
        params_for_star=params_for_star,
        json_out=json_out,
        star_out=star_out,
    )
    # Optionally apply corrections
    if apply_corrections:
        # Decide corrected output path
        corrected_path = Path(output_corrected) if output_corrected else base_path.with_name(f"{data_id}_corrected").with_suffix(".str")
        apply_corrections_from_report(
            report=report,
            input_star=str(Path(str_file)),
            output_star=corrected_path,
            data_id=data_id,
            atoms=correction_atoms,
            release_author=release_author,
        )

    return report


def run_lacs(str_file: str, method: str, data_id: str, rc_model: Optional[Sequence[str] | str] = None,
             outdir: Optional[Path] = None, plots: bool = True, cutoff_k: float = 5.0,
             min_per_side: int = MIN_PER_SIDE_DEFAULT,
             write_format: str = "json",  # {'json','star','both'}
             json_out: Optional[Path] = None,
             star_out: Optional[Path] = None,
             params_for_star: Optional[Dict[str, Any]] = None,  # if None, we'll build one internally
             apply_corrections: bool = False,
             correction_atoms: Sequence[str] = ("CA", "CB", "C", "N"),
             release_author: str = "BMRB",
             output_corrected: Optional[Path] = None) -> Dict[str, Dict]:
    """Run the selected robust method over an NMR-STAR file.

    :param str_file: Path to NMR-STAR file.
    :param method: Robust regression method to use.
    :param data_id: Identifier for the dataset/entry.
    :param rc_model: Random-coil model alias(es), e.g. wis wan; omit for average of all.
    :param outdir: Output directory for plots.``None`` defaults to ``./lacs_output``.
    :param plots: Whether to generate plots.
    :param cutoff_k: Outlier cutoff multiplier.(k in |r|/(k·MAD)).
    :param min_per_side: Minimum number of points required on each sign side.
    :return: Dictionary of results, keyed by list ID.

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
                    fr, draws = _fit_atom_bayes(xvals, yvals, tg, min_per_side=min_per_side)
                    alpha_samples[atom] = draws
                else:
                    fr = _fit_atom_by_method(method, xvals, yvals, tg, min_per_side=min_per_side)
                fits[atom] = fr

        if fits:
            maybe_plot_all(fits, outdir, data_id, method, plots, list_id,cutoff_k=cutoff_k)
            if method == "bayes":
                results[list_id] = collect_and_report_bayes(fits, alpha_samples, cutoff_k=cutoff_k)
            else:
                results[list_id] = collect_and_report(fits, cutoff_k=cutoff_k)
    # ---------- Writing (JSON / STAR / both) ----------
    base_dir = Path.cwd() if outdir is None else Path(outdir)
    base_dir.mkdir(parents=True, exist_ok=True)
    base_path = base_dir / f"{data_id}_{method}"
    report = results
    # If caller didn't pass STAR metadata, build a minimal one from our args
    if params_for_star is None:
        params_for_star = dict(
            str_file=str_file,
            method=method,
            data_id=data_id,
            rc_model=rc_model if rc_model is not None else "",
            outdir=str(base_dir),
            plots=plots,
            cutoff_k=cutoff_k,
            min_per_side=min_per_side,
        )

    written = write_report(
        report=report,
        base_path=base_path,
        write_format=write_format,
        params_for_star=params_for_star,
        json_out=json_out,
        star_out=star_out,
    )
    # Optionally apply corrections
    if apply_corrections:
        # Decide corrected output path
        corrected_path = Path(output_corrected) if output_corrected else base_path.with_name(f"{data_id}_corrected").with_suffix(".str")
        apply_corrections_from_report(
            report=report,
            input_star=str(Path(str_file)),
            output_star=corrected_path,
            data_id=data_id,
            atoms=correction_atoms,
            release_author=release_author,
        )

    return report


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

    :param argv: Command-line arguments (defaults to ``None``).


    Side Effects:

    - Generates optional Plotly plots (HTML and PDF if ``kaleido`` is installed).
    - Writes a JSON report to ``--json-out`` or to a sensible default.

    Examples:

    .. code-block:: bash

        python lacs_unified.py entry.str --method theilsen --data-id demo --out figs --json-out demo_theilsen.json

    """
    import argparse

    p = argparse.ArgumentParser(description="LACS unified CLI for robust linear fits.")
    p.add_argument("star_file", help="Path to NMR-STAR .str file")
    p.add_argument("--method", default='bayes',
                   choices=["tukey","theilsen","ransac","quantile","bayes"],
                   help="Robust regression method to use (default Bayes)")
    p.add_argument("--data-id", default="LACS")
    p.add_argument("--rc-model", nargs="*", default=None,
                   help="Random-coil model alias(es), e.g. wis wan; omit for average of all")
    p.add_argument("--out", type=Path, default=None, help="Output directory for plots")
    p.add_argument("--no-plots", action="store_true", help="Disable plotting entirely")
    p.add_argument('--min-per-side', type=int, default=MIN_PER_SIDE_DEFAULT,
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


    args = p.parse_args(argv)
    rc_model = None if args.rc_model == [] else (args.rc_model if args.rc_model is not None else None)
    report = run_lacs(args.star_file, method=args.method, data_id=args.data_id,
                      rc_model=rc_model, outdir=args.out, plots=not args.no_plots, cutoff_k=args.cutoff_k, min_per_side=args.min_per_side,
                      write_format=args.out_format,
                      json_out=Path(args.json_out) if args.json_out else None,
                      star_out=Path(args.star_out) if args.star_out else None,
                      params_for_star=vars(args),  # capture CLI args as metadata
                      apply_corrections=bool(args.apply_corrections),
                      correction_atoms=args.correction_atoms,
                      release_author=args.release_author,
                      output_corrected=(Path(args.output_corrected) if args.output_corrected else None),
                      )

    # json_path = args.json_out or _default_json_path(args.out, args.data_id, args.method)
    # json_path.parent.mkdir(parents=True, exist_ok=True)
    # with open(json_path, "w", encoding="utf-8") as f:
    #     json.dump(report, f, indent=2)
    # if args.apply_offsets:
    #     if apply_selected_offsets_and_note is None:
    #         raise SystemExit("Cannot apply offsets: apply_lacs_correction.py not importable.")
    #     if not args.output_corrected:
    #         corrected_path = json_path.parent / f'{args.data_id}_corrected,str'
    #         corrected_path.parent.mkdir(parents=True, exist_ok=True)
    #         #raise SystemExit("--output-corrected is required with --apply-offsets.")
    #     else:
    #         corrected_path = args.output_corrected
    #     for list_id in report:
    #         try:
    #             offsets_uc = _extract_offsets_for_list(report, list_id)
    #         except KeyError as e:
    #             # Graceful fallback: try reading the freshly written JSON (identical content)
    #             with open(json_path, "r", encoding="utf-8") as f:
    #                 report2 = json.load(f)
    #             offsets_uc = _extract_offsets_for_list(report2, list_id)
    #
    #         atoms_use = _normalize_atoms(args.atoms)
    #         parts = [f"{a}={offsets_uc[a]:+g}" for a in atoms_use]
    #         details = (
    #                 f"LACS correction applied to list_id {list_id} ({args.data_id}): "
    #                 + (", ".join(parts) if parts else "none")
    #                 + ". Source: LACS (same run)."
    #         )
    #         if len(report)>1:
    #             list_ids = list(report.keys())
    #             if list_ids.index(list_id)>0:
    #                 counts = apply_selected_offsets_and_note(
    #                     input_path=args.star_file,
    #                     output_path=corrected_path,
    #                     list_id=int(list_id),
    #                     offsets=offsets_uc,
    #                     atoms=atoms_use,
    #                     release_author=args.release_author,
    #                     release_details=details,
    #                 )
    #             else:
    #                 counts = apply_selected_offsets_and_note(
    #                     input_path=args.star_file,
    #                     output_path=corrected_path,
    #                     list_id=int(list_id),
    #                     offsets=offsets_uc,
    #                     atoms=atoms_use,
    #                     release_author=args.release_author,
    #                     release_details=details,
    #                 )
    #         else:
    #             counts = apply_selected_offsets_and_note(
    #                 input_path=args.star_file,
    #                 output_path=corrected_path,
    #                 list_id=int(list_id),
    #                 offsets=offsets_uc,
    #                 atoms=atoms_use,
    #                 release_author=args.release_author,
    #                 release_details=details,
    #             )
    #         print("Applied offsets to file:")
    #         for a in atoms_use:
    #             print(f"  {a}: {counts[a]}")
    #         print(f"Total updated: {counts['total']}")
    #     print(f"Wrote corrected file: {corrected_path}")




if __name__ == "__main__":
    main()