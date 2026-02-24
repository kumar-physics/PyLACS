Theory
======

LACS linear relationship
------------------------

For a given nucleus, pylacs fits:

.. math::

   y = b + m x,

where:

- :math:`y` is a nucleus-specific deviation from random coil (:math:`\Delta\delta_A`)
- :math:`x = \Delta\delta_{CA} - \Delta\delta_{CB}`
- :math:`b` and :math:`m` are intercept and slope parameters

Piecewise behavior
------------------

pylacs performs separate fits on each sign-partition of :math:`x`:

- :math:`x \ge 0`
- :math:`x < 0`

This mirrors empirical observations that the relationship differs depending on local structural context.

Robustness and Bayesian uncertainty
-----------------------------------

Outliers may arise from assignment errors, missing corrections, or conformational heterogeneity.
Robust regression methods reduce sensitivity to such points.

When ``--method bayes`` is used, pylacs additionally reports posterior uncertainty for offsets,
derived from posterior samples of per-side intercepts (see :func:`pylacs.lacs.collect_and_report_bayes`).
