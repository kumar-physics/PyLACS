Overview
========

Scientific scope
----------------

pylacs estimates per-atom referencing offsets from an NMR-STAR chemical shift list by:

1. computing deviations from random-coil references, :math:`\Delta\delta`,
2. forming a predictor from carbon shifts,
3. fitting robust linear models per nucleus and per sign-partition,
4. reporting offsets and outlier probabilities.

Atoms analyzed
--------------

The analysis is performed for the following atoms when present:

- ``CA``
- ``CB``
- ``C``
- ``N``
- ``HA``

Random-coil reference models
----------------------------

Random-coil values are provided by :class:`pylacs.random_coil.RandomCoil`.
The following model families are available (with accepted aliases):

- Wishart (``wis`` / ``wishart``)
- Wang (``wan`` / ``wang``)
- Lukhin (``luk`` / ``lukhin``)
- Schwarzinger (``sch`` / ``schwarzinger``)
- Poulsen (temperature-dependent; ``pou`` / ``poulsen``)

If no model is specified, :meth:`pylacs.random_coil.RandomCoil.get_value` averages across **all** models.

Primary entry points
--------------------

- CLI: :mod:`pylacs.cli` (installed script ``pylacs``)
- Library API: :func:`pylacs.lacs.run_lacs`
- Correction utility: :func:`pylacs.apply_lacs_correction.apply_selected_offsets_and_note`
