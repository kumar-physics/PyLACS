pylacs
======

**PyLACS** is a Python implementation of **Linear Analysis of Chemical Shifts (LACS)** :cite:`Wang2005` for identifying and
quantifying **systematic chemical-shift referencing offsets** in biomolecular NMR datasets stored as
NMR-STAR (``.str``) :cite:`Ulrich2019` chemical shift files.

pylacs provides:

- a command-line program ``pylacs`` (recommended entry point),
- a programmatic API centered around :func:`pylacs.lacs.run_lacs`,
- an optional post-processing utility that applies estimated offsets back to an NMR-STAR file
  and appends a Release record (:func:`pylacs.apply_lacs_correction.apply_selected_offsets_and_note`).

.. note::
   This documentation corresponds to version |version|.

.. toctree::
   :maxdepth: 2
   :caption: User Guide

   overview
   getting_started
   usage
   methodology
   theory
   outputs
   examples
   validation
   release-notes
   citation

.. toctree::
   :maxdepth: 2
   :caption: API Reference

   api/modules
