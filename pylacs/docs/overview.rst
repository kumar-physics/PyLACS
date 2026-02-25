Overview
========

Scientific scope
----------------

pylacs estimates per-atom type referencing offsets from an NMR-STAR chemical shift list by:

1. computing deviations from random-coil references, :math:`\Delta\delta`,
2. forming a predictor from carbon shifts,
3. fitting robust linear models per nucleus and per sign-partition,
4. reporting offsets and outlier probabilities.

Requirement
____________
The chemical shift data should contain ``CA`` and ``CB`` chemical shifts, otherwise it is not possible
to calculate the reference independent quantity :math:`\Delta \delta CA - \Delta \delta CB` where
:math:`\Delta \delta A = \delta_{Observed} A-\delta_{Random Coil} A`

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

- Wishart (``wis`` / ``wishart``) :cite:`Wishart1995`
- Wang (``wan`` / ``wang``) :cite:`Wang2002`
- Lukhin (``luk`` / ``lukhin``) :cite:`Lukin1997`
- Schwarzinger (``sch`` / ``schwarzinger``) :cite:`Schwarzinger2000`
- Poulsen (temperature-dependent; ``pou`` / ``poulsen``) :cite:`Kjaergaard2011`

If no model is specified, :meth:`pylacs.random_coil.RandomCoil.get_value` averages across **all** models.

Primary entry points
--------------------

- CLI: :mod:`pylacs.cli` (installed script ``pylacs``)
- Library API: :func:`pylacs.lacs.run_lacs`
- Correction utility: :func:`pylacs.apply_lacs_correction.apply_selected_offsets_and_note`
