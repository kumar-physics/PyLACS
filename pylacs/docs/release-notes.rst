Release Notes
=============

Version 1.0.1
-------------

- Requirements updated to numpy version 2.x
- Unused CLI option removed
- Documentation updated

Version 1.0.0
-------------

- Unified implementation of LACS offset estimation in :mod:`pylacs.lacs`.
- Five regression methods: Tukey RLM, Theil–Sen, RANSAC, Quantile regression, Bayesian Student-t.
- Structured JSON reporting with per-residue outlier probabilities and (for Bayes) posterior summaries.
- Optional application of offsets to an NMR-STAR file with appended Release record
  (:mod:`pylacs.apply_lacs_correction`).
