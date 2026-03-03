PyLACS
======

**PyLACS** is a Python implementation of **Linear Analysis of Chemical Shifts (LACS)** :cite:`Wang2005` for identifying and
quantifying **systematic chemical-shift referencing offsets** in biomolecular NMR datasets stored as
NMR-STAR (``.str``) :cite:`Ulrich2019` chemical shift files. PyLAC is developed and maintained by **Biological Magnetic
Resonance data Bank (BMRB)** :cite:`Ulrich2007,Romero2020,Hoch2023`

PyLACS provides:

- a command-line program ``pylacs`` (recommended entry point),
- a programmatic API centered around :func:`pylacs.lacs.run_lacs`,
- an option to chose different robust linear fit algorithm
- an option to chose different Random Coil chemical shifts or to use average value from more than one Random Coil Shifts
- an optional post-processing utility that applies estimated offsets back to an NMR-STAR file
  and appends a Release record (:func:`pylacs.apply_lacs_correction.apply_selected_offsets_and_note`).

.. note::
   This documentation corresponds to version |version|.

.. toctree::
   :maxdepth: 2
   :caption: User Guide

   overview
   getting_started
   theory
   methodology
   usage
   outputs
   examples
   validation
   release-notes
   citation

.. toctree::
   :maxdepth: 2
   :caption: API Reference

   api
