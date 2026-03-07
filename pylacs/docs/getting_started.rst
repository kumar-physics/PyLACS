Quick Start
===========

This section provides a minimal, reproducible workflow for running LACS :cite:`Wang2005`  analysis with pylacs.
All examples below assume an NMR-STAR chemical shift file (``.str``).

Installation
------------

Install from PyPI:

.. code-block:: bash

   pip install pylacs

Local development install (editable):

.. code-block:: bash

   git clone <YOUR_REPO_URL>
   cd pylacs
   pip install -e .

Minimal CLI run
---------------

Run the analysis using the default method (``bayes``) and default metadata id (``LACS``):

.. code-block:: bash

   pylacs path/to/ENTRY.str

Specify a dataset identifier for filenames and metadata:

.. code-block:: bash

   pylacs path/to/ENTRY.str --data-id 2KOC

Choose a robust regression method:

.. code-block:: bash

   pylacs ENTRY.str --method tukey
   pylacs ENTRY.str --method theilsen
   pylacs ENTRY.str --method ransac
   pylacs ENTRY.str --method quantile
   pylacs ENTRY.str --method bayes

Select random-coil reference model(s)
-------------------------------------

By default, random-coil values are averaged across all models in
:class:`pylacs.random_coil.RandomCoil`.

To restrict the reference to specific models, use ``--rc-model`` with one or more aliases:

.. code-block:: bash

   pylacs ENTRY.str --rc-model wis
   pylacs ENTRY.str --rc-model wis wan
   pylacs ENTRY.str --rc-model kja



Where outputs go
----------------

- JSON/STAR reports are written to the current working directory by default,
  with names like ``<data-id>_<method>.json`` and ``<data-id>_<method>.str``.
- Diagnostic plots are written to ``./lacs_output`` by default (unless ``--out`` is provided).

See :doc:`outputs` for details.
