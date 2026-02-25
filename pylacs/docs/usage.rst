Usage
=====

Command-line interface
----------------------

The installed command is:

.. code-block:: bash

   pylacs STAR_FILE.str [options]

The CLI is implemented in :mod:`pylacs.cli` and calls :func:`pylacs.lacs.run_lacs`.

Core options
^^^^^^^^^^^^

``--method {tukey,theilsen,ransac,quantile,bayes}``
  Robust regression method. Default: ``bayes``.

``--data-id TEXT``
  Identifier propagated into filenames and output metadata. Default: ``LACS``.

``--rc-model ALIAS [ALIAS ...]``
  One or more random-coil model aliases (e.g., ``wis``, ``wan``, ``luk``, ``sch``, ``pou``).
  If omitted, random-coil values are averaged across all models.

``--min-per-side INT``
  Minimum number of points required *on each sign side* of :math:`x` to fit both sides.
  Default: 5.

``--cutoff-k FLOAT``
  Multiplier :math:`k` used in outlier scoring from robust-scaled residuals. Default: 5.0.

Plotting and output control
^^^^^^^^^^^^^^^^^^^^^^^^^^^

``--out PATH``
  Directory for plots (and used as the report base path if explicit ``--json-out``/``--star-out`` are not set).

``--no-plots``
  Disable plotting.

``--out-format {json,star,both}``
  Write JSON, STAR, or both. Default: ``both``.

``--json-out PATH``
  Explicit JSON output path.

``--star-out PATH``
  Explicit STAR output path.

Applying corrections
^^^^^^^^^^^^^^^^^^^^

pylacs can apply the estimated offsets to the input NMR-STAR file and write a corrected file:

``--apply-corrections``
  Apply offsets to selected atoms in the chemical shift loop(s) and append a Release record.

``--correction-atoms ATOM [ATOM ...]``
  Atoms to correct when applying offsets. Default: ``CA CB C N``.

``--release-author TEXT``
  Recorded author in the Release record. Default: ``BMRB``.

``--output-corrected PATH``
  Output filename for the corrected STAR file.

