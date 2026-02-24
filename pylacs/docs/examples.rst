Examples
========

Example 1: Default Bayesian analysis
------------------------------------

.. code-block:: bash

   pylacs ENTRY.str --data-id ENTRY

Example 2: Tukey robust regression
----------------------------------

.. code-block:: bash

   pylacs ENTRY.str --method tukey --data-id ENTRY --out plots/

Example 3: Restrict random-coil reference model
-----------------------------------------------

.. code-block:: bash

   pylacs ENTRY.str --rc-model wis wan --data-id ENTRY

Example 4: Apply offset corrections to a STAR file
--------------------------------------------------

.. code-block:: bash

   pylacs ENTRY.str --apply-corrections --correction-atoms CA CB C N \
     --release-author "BMRB" --output-corrected ENTRY_corrected.str
