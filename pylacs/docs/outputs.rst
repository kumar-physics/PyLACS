Outputs
=======

JSON report structure
---------------------

The output of :func:`pylacs.lacs.run_lacs` is a dictionary keyed by chemical-shift list ID.
Each list contains:

- ``offsets``: per-atom offsets (rounded to 3 decimals)
- ``offsets_split``: per-atom offsets for ``pos`` and ``neg`` sign partitions
- ``outliers``: per-atom list of residue records, each with:

  - ``residue_key``: entity, assembly, index, comp
  - ``x`` and ``y`` (the predictor and response used in the fit)
  - ``residual``
  - ``prob`` (soft outlier probability)
  - ``flag`` (binary indicator based on the cutoff)

- ``meta``: includes ``cutoff_k`` and ``used_side`` per atom

For the Bayesian method, additional keys are included:

- ``offsets_bayes``: mean, 95% credible interval, and SD of the offset
- ``offsets_bayes_sides``: posterior summaries per side (pos/neg)

See :func:`pylacs.lacs.collect_and_report` and :func:`pylacs.lacs.collect_and_report_bayes`.

File naming
-----------

If explicit output paths are not provided:

- JSON: ``<data-id>_<method>.json``
- STAR: ``<data-id>_<method>.str``

These defaults are computed relative to the current working directory or ``--out``.

Plots
-----

If plotting is enabled and Plotly is available, pylacs writes per-atom diagnostic plots.
If ``--out`` is not provided, plot files are written to ``./lacs_output`` (created if needed).

See :func:`pylacs.lacs.plot_all`.
