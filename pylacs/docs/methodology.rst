Methodology
===========

Input parsing (NMR-STAR)
------------------------

The input NMR-STAR file is parsed using **pynmrstar**.
Chemical shifts are extracted per chemical-shift list and assembled into a per-residue mapping.

See:

- :func:`pylacs.lacs.read_star`
- :func:`pylacs.lacs._star_tag_value`
- :func:`pylacs.lacs._star_loop_cell`

Random-coil referencing and :math:`\Delta\delta`
------------------------------------------------

For each residue and atom, pylacs computes deviations from random-coil references:

.. math::

   \Delta\delta_A = \delta_A^{obs} - \delta_A^{RC}.

Random-coil values are obtained from :class:`pylacs.random_coil.RandomCoil` via
:meth:`pylacs.random_coil.RandomCoil.get_value`.

Predictor definition
--------------------

The LACS predictor is constructed when both ``CA`` and ``CB`` are present for a residue:

.. math::

   x = \Delta\delta_{CA} - \Delta\delta_{CB}.

For nuclei ``C``, ``N``, and ``HA``, pylacs stores the corresponding :math:`x` values alongside each
:math:`y = \Delta\delta_A` value.

Implementation details:

- :func:`pylacs.lacs.compute_deltas`

Sign-partitioned regression
---------------------------

Because empirical behavior differs by the sign of :math:`x`, pylacs partitions data into:

- positive side: :math:`x \ge 0`
- negative side: :math:`x < 0`

Each side is fit independently (if enough points are present).

Offsets and sign convention
---------------------------

For each atom, pylacs records per-side intercept-based offsets as:

- ``offsets_split[atom]["pos"] = -intercept_pos``
- ``offsets_split[atom]["neg"] = -intercept_neg``

The final offset is:

- if both sides exist: :math:`-\frac{1}{2}(b_{pos} + b_{neg})`
- if only one side exists: :math:`-b_{side}`
- if neither side exists: 0.0

See :func:`pylacs.lacs.collect_and_report`.

Robust regression methods
-------------------------

Per-side fits are provided by:

- Tukey biweight RLM (statsmodels): :func:`pylacs.lacs.fit_side_rlm_tukey`
- Theil–Sen (scikit-learn): :func:`pylacs.lacs.fit_side_theilsen`
- RANSAC (scikit-learn): :func:`pylacs.lacs.fit_side_ransac`
- Quantile regression at :math:`\tau=0.5` (statsmodels): :func:`pylacs.lacs.fit_side_quantile`
- Bayesian Student-t regression (PyMC): :func:`pylacs.lacs.fit_side_bayes_t`

Outlier probabilities
---------------------

Residuals are robustly scaled using MAD (median absolute deviation) and mapped to a smooth probability
using a logistic function.

See:

- :func:`pylacs.lacs.mad`
- :func:`pylacs.lacs.logistic_prob`
- :func:`pylacs.lacs.outlier_stats`
