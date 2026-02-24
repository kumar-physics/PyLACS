Validation and Quality Control
==============================

Minimum data requirements
-------------------------

pylacs requires residues with both ``CA`` and ``CB`` assignments in order to compute

.. math::

   x = \Delta\delta_{CA} - \Delta\delta_{CB}.

For two-sided fitting, pylacs additionally requires at least ``--min-per-side`` points on each sign side.

Interpreting outliers
---------------------

Outlier probabilities are soft scores derived from robust-scaled residuals and should be interpreted as
prioritization for manual inspection rather than definitive classification.

Recommended QC steps
--------------------

- Check the number of points per side for each nucleus.
- Compare offsets across methods (e.g., ``tukey`` vs ``theilsen`` vs ``bayes``).
- Inspect high-probability outliers to identify systematic issues in assignments or referencing.
