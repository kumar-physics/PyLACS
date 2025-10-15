# PyLACS : Python version Linear Analysis of Chemical Shifts(LACS)



[![PyPI version](https://badge.fury.io/py/pylacs.svg)](https://pypi.org/project/pylacs/)
[![Python Versions](https://img.shields.io/pypi/pyversions/pylacs.svg)](https://pypi.org/project/pylacs/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

[//]: # ([![Build Status]&#40;https://github.com/<your-org>/pylacs/actions/workflows/tests.yml/badge.svg&#41;]&#40;https://github.com/<your-org>/pylacs/actions&#41;)

[//]: # ([![DOI]&#40;https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg&#41;]&#40;https://doi.org/10.5281/zenodo.XXXXXXX&#41;)

**PyLACS** is a Python package for **detecting and correcting chemical shift referencing errors** in protein NMR spectroscopy datasets, using robust linear regression methods.  
It provides both a **Python API** and a **command-line interface (CLI)** for use in research pipelines, HPC batch jobs, and interactive exploration.

---

## ✨ Features

- Detect and correct chemical shift referencing errors for **Cα, Cβ, C, N, and H** nuclei  
- Robust regression methods:
  - Tukey biweight  
  - Theil–Sen estimator  
  - RANSAC  
  - Quantile regression  
  - Bayesian regression (via PyMC + ArviZ)  
- Outlier detection and offset estimation  
- Optional **plot generation** (scatter, regression fits)  
- Works directly with **BMRB NMR-STAR (.str) files**  
- Provides both a **CLI** (`pylacs`) and a **Python API**  

---

## 📦 Installation

Stable release from PyPI:

```bash
pip install pylacs
```

---
## 🚀 Quickstart
```bash
pylacs myfile.str --method tukey --out results/
````

This will compute offsets, generate validation reports, and (if plotting is enabled) save plots in results folder

For help

```bash
pylacs --help
```
## 📖 Citation

f you use PyLACS in your research, please cite:

Wang, L., Eghbalnia, H. R., & Markley, J. L. (2005).
Probabilistic approach to determining protein backbone torsion angles from NMR chemical shifts.
Journal of Biomolecular NMR, 32(1), 13–22.

Baskaran, K., et al. (2025).
PyLACS: Python-based Linear Analysis of Chemical Shifts pipeline.
(manuscript in preparation)
