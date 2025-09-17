#!/usr/bin/env python3
"""
Random coil chemical-shift utilities.

This module provides a :class:`RandomCoil` class that exposes several published
random-coil chemical shift reference tables and convenient helpers to access
per-residue/atom values, compute model averages, and validate inputs.

Key features
------------
- Access four literature models ("wishart", "wang", "lukhin", "schwarzinger").
- Query a specific residue/atom value with one- or three-letter residue codes.
- Compute averages across an arbitrary subset of models or across **all** models.
- Strict input validation with clear error messages.

Examples
--------
>>> rc = RandomCoil()
>>> rc.get_value('HIS', 'C', ['wis', 'wan'])
176.02  # example output

Design notes
------------
- All API surface is typed to improve IDE support and static analysis.
- Public API is intentionally small: prefer :meth:`RandomCoil.get_value` and
  :meth:`RandomCoil.get_average` over accessing model dicts directly.

"""

from __future__ import annotations

from typing import Dict, Iterable, List, Mapping, Sequence, Tuple, Union
import math
import numpy as np
import warnings


ResidueCode = str
AtomName = str
ModelName = str
Scalar = float


class RandomCoil:
    """Access random-coil chemical shift reference values.

    Parameters
    ----------
    None

    Notes
    -----
    - Supported atom names (case-insensitive): ``N, C, CA, CB, H, HA``.
    - For glycine, ``HA2``/``HA3`` are normalized to ``HA``.
    - Model name aliases (short/long) are accepted, e.g. ``'wis'`` or
      ``'wishart'``.

    """

    # --- residue code maps -------------------------------------------------
    _THREE_TO_ONE: Mapping[str, str] = {
        'ILE': 'I', 'GLN': 'Q', 'GLY': 'G', 'GLU': 'E', 'CYS': 'C',
        'ASP': 'D', 'SER': 'S', 'LYS': 'K', 'PRO': 'P', 'ASN': 'N',
        'VAL': 'V', 'THR': 'T', 'HIS': 'H', 'TRP': 'W', 'PHE': 'F',
        'ALA': 'A', 'MET': 'M', 'LEU': 'L', 'ARG': 'R', 'TYR': 'Y',
    }
    _ONE_TO_THREE: Mapping[str, str] = {v: k for k, v in _THREE_TO_ONE.items()}

    # --- atom order shared by all model tables ----------------------------
    _ATOMS: Tuple[str, ...] = ('N', 'C', 'CA', 'CB', 'H', 'HA')

    # ------------------------------- API -----------------------------------
    def atoms(self) -> List[str]:
        """Return the list of supported atom names in canonical order."""
        return list(self._ATOMS)

    def get_names(self, short_or_long: str = 'short') -> List[str]:
        """Return the model names in either short or long form.

        Parameters
        ----------
        short_or_long:
            ``'short'`` -> ``['wis','wan','luk','sch','pou']``;
            ``'long'``  -> ``['wishart','wang','lukhin','schwarzinger','poulsen']``

        """
        if short_or_long == 'short':
            return ['wis', 'wan', 'luk', 'sch', 'pou']
        if short_or_long == 'long':
            return ['wishart', 'wang', 'lukhin', 'schwarzinger', 'poulsen']
        raise ValueError("short_or_long must be 'short' or 'long'")

    # ------------------------------ helpers --------------------------------
    def _normalize_residue(self, res: ResidueCode) -> str:
        res = res.upper()
        if len(res) == 3:
            if res not in self._THREE_TO_ONE:
                raise ValueError(f"Unknown three-letter code: {res}")
            return self._THREE_TO_ONE[res]
        if len(res) == 1:
            if res not in self._ONE_TO_THREE:
                raise ValueError(f"Unknown one-letter code: {res}")
            return res
        raise ValueError(f"Residue code must be 1 or 3 letters, got: {res}")

    def _normalize_atom(self, atom: AtomName, res_one: str) -> str:
        atom = atom.upper()
        if res_one == 'G' and atom in {'HA2', 'HA3'}:
            atom = 'HA'
        if atom not in self._ATOMS:
            raise ValueError(f"Unknown atom '{atom}'. Supported: {', '.join(self._ATOMS)}")
        return atom

    def _select_models(self, rc_name: Union[None, str, Sequence[str]]) -> Dict[str, Mapping[str, Sequence[Scalar]]]:
        """Return a mapping of model_name -> table for the requested models.

        If ``rc_name`` is ``None``, all models are used (for averaging). If a
        string is provided, it can be any accepted alias. If a sequence is
        provided, each may be an alias.
        """
        alias_map = {
            'wis': 'wishart', 'wishart': 'wishart',
            'wan': 'wang', 'wang': 'wang',
            'luk': 'lukhin', 'lukhin': 'lukhin',
            'sch': 'schwarzinger', 'schwarzinger': 'schwarzinger',
            'pou': 'poulsen', 'poulsen': 'poulsen'
        }
        if rc_name is None:
            names = ['wishart', 'wang', 'lukhin', 'schwarzinger','poulsen']
        elif isinstance(rc_name, str):
            if rc_name.lower() not in alias_map:
                raise ValueError(f"Unknown random-coil model: {rc_name}")
            names = [alias_map[rc_name.lower()]]
        else:
            names = []
            for nm in rc_name:
                key = nm.lower()
                if key not in alias_map:
                    raise ValueError(f"Unknown random-coil model: {nm}")
                names.append(alias_map[key])
        models = {
            'wishart': self.wishart,
            'wang': self.wang,
            'lukhin': self.lukhin,
            'schwarzinger': self.schwarzinger,
            'poulsen': self.poulsen
        }
        return {n: models[n] for n in names}

    # ------------------------------- public ---------------------------------
    def get_value(self, res: ResidueCode, atom: AtomName,
                  rc_name: Union[None, str, Sequence[str]] = None,temp: float=25.0) -> Scalar:
        """Return the random-coil value (ppm) for a residue/atom.

        Parameters
        ----------
        res:
            Residue code (one- or three-letter). Case-insensitive.
        atom:
            Atom name (e.g., ``'CA'``, ``'N'``). Case-insensitive.
        rc_name:
            - ``None``: average across **all** models
            - ``str``: a single model (alias accepted, e.g., ``'wis'``)
            - ``Sequence[str]``: average across the specified models
        temp:
            Temperature in C It will be used only when the rc_name is pou or poulsen

        Returns
        -------
        float
            The random-coil chemical shift value (ppm).

        Raises
        ------
        ValueError
            If the residue code, atom name, or model name is invalid.
        """
        r1 = self._normalize_residue(res)
        atom = self._normalize_atom(atom, r1)

        # Select models
        selected = self._select_models(rc_name)

        # Collect values across selected models
        idx = self._ATOMS.index(atom)
        values: List[float] = []
        for _, table in selected.items():
            val = table[r1][idx]
            values.append(math.nan if val is None else float(val))


        # Average across models (handles NaNs)
        arr = np.array(values, dtype=float)
        if rc_name == 'pou' or rc_name == ['pou'] or rc_name == 'poulsen' or rc_name == ['poulsen']:
            coeff: List[float] = []
            for _, table in {'ce': self.poulsen_temp_coeff}.items():
                val = table[r1][idx]
                coeff.append(math.nan if val is None else float(val))
                coef = np.array(coeff, dtype=float)
                cs_t = arr + coef * (temp - 25.0)
                return round(float(np.nanmean(cs_t)),2)
        else:
            return float(np.nanmean(arr))


    def get_average(self, name_list: Sequence[str] | None = None) -> Mapping[str, List[float]]:
        """
        Return per-residue averages across the given models (default: all).

        - If a residue is missing from a model's table, it is treated as all-NaN for that model.
        - Averages are computed with np.nanmean per atom, so missing entries are ignored.
        """
        selected = self._select_models(None if name_list is None else name_list)

        # All residues that appear in at least one selected model
        residues = set().union(*(tbl.keys() for tbl in selected.values()))

        n_atoms = len(self.atoms())  # expected vector length, e.g. 6: N, C, CA, CB, H, HA
        avg: Dict[str, List[float]] = {}

        for r in sorted(residues):
            mats = []
            for tbl in selected.values():
                vec = tbl.get(r)
                if vec is None:
                    mats.append([math.nan] * n_atoms)
                else:
                    mats.append([math.nan if v is None else float(v) for v in vec])

            arr = np.asarray(mats, dtype=float)  # shape: (n_models, n_atoms)

            # Optional: suppress RuntimeWarning “Mean of empty slice” when a whole column is NaN
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
                colmeans = np.nanmean(arr, axis=0)

            avg[r] = [float(round(x, 2)) if not np.isnan(x) else math.nan for x in colmeans]

        return avg

    # ------------------------------ data tables -----------------------------
    # The following tables are copied from the user's original module.
    # Atom order: N, C, CA, CB, H, HA
    wishart: Mapping[str, Sequence[Scalar]] = {
        "I": [121.7, 176.4, 61.1, 38.8, 8.00, 4.17],
        "V": [120.5, 176.3, 62.2, 32.9, 8.03, 4.12],
        "D": [121.4, 176.3, 54.2, 41.1, 8.34, 4.64],
        "N": [119.0, 175.2, 53.1, 38.9, 8.40, 4.74],
        "F": [120.9, 175.8, 57.7, 39.6, 8.30, 4.62],
        "H": [118.2, 174.1, 55.0, 29.0, 8.42, 4.73],
        "W": [122.2, 176.1, 57.5, 29.60, 8.25, 4.66],
        "Y": [120.8, 175.9, 57.9, 38.8, 8.12, 4.55],
        "K": [121.6, 176.6, 56.2, 33.1, 8.29, 4.32],
        "L": [122.6, 177.6, 55.1, 42.4, 8.16, 4.34],
        "M": [120.7, 176.3, 55.4, 32.9, 8.28, 4.48],
        "Q": [120.6, 176.0, 55.7, 29.4, 8.32, 4.34],
        "R": [121.3, 176.3, 56.0, 30.9, 8.23, 4.34],
        "E": [121.7, 176.6, 56.6, 29.9, 8.42, 4.35],
        "T": [116.0, 174.7, 61.8, 69.8, 8.15, 4.35],
        "C": [118.9, 174.6, 58.2, 28.0, 8.43, 4.71],
        "S": [116.6, 174.6, 58.3, 63.8, 8.31, 4.47],
        "A": [125.0, 177.8, 52.5, 19.1, 8.24, 4.32],
        "G": [109.1, 174.9, 45.1, None, 8.33, 3.96],
        "P": [None, 177.3, 63.3, 32.1, None, 4.42],
        "B": [118.6, 174.6, 55.4, 41.1, 8.43, 4.71],
    }

    wang: Mapping[str, Sequence[Scalar]] = {
        "I": [120.58, 175.52, 60.79, 38.43, 7.94, 4.18],
        "V": [119.91, 175.66, 62.00, 32.35, 7.98, 4.13],
        "D": [120.37, 176.00, 54.00, 40.78, 8.31, 4.62],
        "N": [118.50, 174.84, 53.00, 38.43, 8.35, 4.66],
        "F": [119.72, 175.46, 57.46, 39.41, 8.09, 4.59],
        "H": [118.92, 174.78, 55.74, 29.50, 8.18, 4.60],
        "W": [120.99, 175.87, 57.54, 29.60, 7.97, 4.60],
        "Y": [119.37, 175.29, 57.64, 38.78, 7.99, 4.56],
        "K": [121.10, 176.15, 56.29, 32.53, 8.17, 4.28],
        "L": [121.57, 176.70, 54.77, 42.14, 8.06, 4.36],
        "M": [120.14, 175.94, 55.43, 32.92, 8.22, 4.47],
        "Q": [119.82, 175.75, 55.89, 29.01, 8.20, 4.29],
        "R": [120.75, 176.01, 56.18, 30.36, 8.21, 4.26],
        "E": [120.62, 176.32, 56.66, 29.87, 8.36, 4.28],
        "T": [113.88, 174.78, 61.30, 68.92, 8.16, 4.44],
        "C": [118.10, 175.11, 58.24, 29.54, 8.10, 4.59],
        "S": [116.00, 174.41, 58.20, 63.75, 8.22, 4.45],
        "A": [123.82, 177.28, 52.46, 18.98, 8.09, 4.31],
        "G": [109.48, 174.01, 45.28, None, 8.37, 3.97],
        "P": [None, 176.62, 63.24, 31.81, None, 4.41],
        "B": [118.7, 175.5, 55.6, 41.2, 8.54, 4.76],
    }

    lukhin: Mapping[str, Sequence[Scalar]] = {
        "G": [110.03, 173.96, 45.41, None, 8.37, 3.97],
        "B": [118.7, 175.5, 55.6, 41.2, 8.54, 4.76],
        "A": [125.10, 177.30, 52.42, 19.03, 8.09, 4.31],
        "S": [116.43, 174.41, 58.27, 64.14, 8.22, 4.45],
        "C": [119.19, 174.84, 58.01, 28.20, 8.43, 4.71],
        "M": [120.39, 175.45, 55.34, 33.00, 8.22, 4.47],
        "K": [121.82, 176.39, 56.59, 32.62, 8.21, 4.26],
        "V": [120.93, 175.79, 62.13, 32.65, 7.97, 4.60],
        "T": [114.17, 174.75, 61.62, 69.83, 8.16, 4.44],
        "I": [121.19, 175.69, 60.98, 38.87, 7.94, 4.18],
        "L": [122.22, 177.15, 54.82, 42.82, 8.06, 4.36],
        "D": [120.51, 176.45, 54.12, 40.83, 8.31, 4.62],
        "N": [119.17, 174.65, 53.22, 38.74, 8.35, 4.66],
        "E": [121.44, 176.27, 56.66, 30.13, 8.36, 4.28],
        "Q": [120.34, 175.54, 55.78, 29.34, 8.20, 4.29],
        "R": [122.42, 176.05, 56.25, 30.56, 8.21, 4.26],
        "H": [120.21, 174.54, 55.78, 29.78, 8.18, 4.60],
        "F": [119.84, 174.79, 57.91, 39.34, 8.09, 4.59],
        "Y": [120.41, 175.80, 57.77, 38.88, 7.99, 4.56],
        "W": [120.19, 175.85, 57.50, 29.09, 7.97, 4.60],
        "P": [136.73, 176.60, 63.27, 32.09, None, 4.41],
    }

    schwarzinger: Mapping[str, Sequence[Scalar]] = {
        "A": [125.0, 178.5, 52.8, 19.3, 8.35, 4.35],
        "B": [118.7, 175.5, 55.6, 41.2, 8.54, 4.76],
        "C": [118.8, 175.3, 58.6, 28.3, 8.44, 4.59],
        "D": [119.1, 175.9, 53.0, 38.3, 8.56, 4.82],
        "E": [120.2, 176.8, 56.1, 29.9, 8.40, 4.42],
        "F": [120.7, 176.6, 58.1, 39.8, 8.31, 4.65],
        "G": [107.5, 174.9, 45.4, None, 8.41, 4.02],
        "H": [118.1, 175.1, 55.4, 29.1, 8.56, 4.79],
        "I": [120.4, 177.1, 61.6, 38.9, 8.17, 4.21],
        "K": [121.6, 177.4, 56.7, 33.2, 8.36, 4.36],
        "L": [122.4, 178.2, 55.5, 42.5, 8.28, 4.38],
        "M": [120.3, 177.1, 55.8, 32.9, 8.42, 4.52],
        "N": [119.0, 176.1, 53.3, 39.1, 8.51, 4.79],
        "P": [None, 177.8, 63.7, 32.2, None, 4.45],
        "Z": [None, None, 63.0, 34.8, None, 4.60],
        "Q": [120.5, 176.8, 56.2, 29.5, 8.44, 4.38],
        "R": [121.2, 177.1, 56.5, 30.9, 8.39, 4.38],
        "S": [115.5, 175.4, 58.7, 64.1, 8.43, 4.51],
        "T": [112.0, 175.6, 62.0, 70.0, 8.25, 4.43],
        "V": [119.3, 177.0, 62.6, 31.8, 8.16, 4.16],
        "W": [122.1, 177.1, 57.6, 29.8, 8.22, 4.70],
        "Y": [120.9, 176.7, 58.3, 38.9, 8.26, 4.58],
    }
    # pH 6.5 1M urea and temperature 25C DOI 10.1007/s10858-011-9472-x
    poulsen: Mapping[str, Sequence[Scalar]] = {
        "A": [124.14, 178.64, 52.77, 19.18, 8.35, 4.37],#
        "C": [118.74, 175.49, 58.64, 28.12, 8.42, 4.59],#
        "D": [120.43, 177.25, 54.42, 41.29, 8.39, 4.65],#
        "E": [121.03, 177.52, 57.02, 30.09, 8.55, 4.32],#
        "F": [120.35, 176.76, 58.06, 39.55, 8.26, 4.65],#
        "G": [109.04, 175.08, 45.37, None,  8.42, 4.02],#
        "H": [118.78, 175.75, 55.88, 29.83, 8.41, 4.73],#
        "I": [120.08, 177.29, 61.57, 38.71, 8.13, 4.21],#
        "K": [120.99, 177.48, 56.54, 32.98, 8.33, 4.37],#
        "L": [121.83, 178.39, 55.48, 42.32, 8.27, 4.38],#
        "M": [119.95, 177.18, 55.71, 32.75, 8.40, 4.54],#
        "N": [118.88, 176.14, 53.63, 38.97, 8.49, 4.78],#
        "P": [None,   178.08, 63.71, 32.09, None, 4.46],#
        "Q": [120.04, 176.88, 56.14, 29.38, 8.43, 4.39],#
        "R": [120.72, 177.23, 56.38, 30.71, 8.35, 4.38],#
        "S": [116.04, 175.43, 58.60, 63.90, 8.41, 4.51],#
        "T": [113.05, 175.64, 62.06, 69.84, 8.24, 4.41], #
        "V": [119.30, 177.21, 62.64, 32.64, 8.13, 4.17],#
        "W": [121.26, 177.18, 57.59, 29.51, 8.16, 4.70],#
        "Y": [120.48, 176.82, 58.24, 38.73, 8.21, 4.58],#
    }
    poulsen_temp_coeff: Mapping[str,Sequence[Scalar]]={
        "A":[-5.3,-7.1,-2.2, 4.7,-9.0, 0.7],
        "C":[-8.2,-2.6,-0.9, 1.3,-7.0, 0.0],
        "D":[-3.9,-4.8,2.8, 6.5,-6.2,-0.1],
        "E":[-3.7,-4.9,0.9, 4.6,-6.5, 0.3],
        "F":[-11.2,-6.9,-4.7, 2.4,-7.5, 0.4],
        "G":[-6.2,-3.2,3.3,None, -9.1, 0.0],
        "H":[3.3, 3.1,7.8, 15.5,  -7.8, 0.4],
        "I":[-12.7,-8.7,-2.0, 4.6,-7.8, 0.4],
        "K":[-7.6,-7.1,-0.8, 2.4,-7.5, 0.4],
        "L":[-2.9, -8.2,1.7 ,4.9,-7.5, 0.1],
        "M":[-6.2,-8.2, 4.1, 9.4,-7.1,-0.5],
        "N":[-3.3, -6.1,2.8, 5.1,-7.0,-2.9],
        "P":[None,-4.0,1.1,-0.2,  None, 0.0],
        "Q":[-6.5,-5.7,2.3, 3.6,-7.2, 0.3],
        "R":[-5.3,-6.9,-1.4, 3.5,-7.1, 0.4],
        "S":[-3.8,-4.7,-1.7, 4.4,-7.6, 0.1],
        "T":[-6.7,-5.2,0.0, 2.2,-7.3, 0.0],
        "V":[-14.2,-8.1,-2.8, 2.5,-7.6, 0.5],
        "W":[-10.1,-7.9,-2.7, 3.1,-7.8, 0.4],
        "Y":[-12.0,-7.7,-5.0, 2.9,-7.7, 0.5]
    }


def _demo() -> None:
    rc = RandomCoil()
    print(rc.get_value('HIS', 'C',['wis','wan']))
    print(rc.get_value('HIS', 'C', 'pou',temp=10))

if __name__ == '__main__':
    _demo()
