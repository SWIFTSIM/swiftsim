.. Special Behavior
   Loic Hausammann, 2020

.. _special_behavior:

Special Behavior
================

Generating new unique IDs
-------------------------

When generating new particles (not transforming them), the code needs to provide new unique IDs.
This is implemented in the file ``space_unique_id.c`` and can be turn on/off in the star formation file ``star_formation_struct.h``
with the variable ``star_formation_need_unique_id``.

The generation of new IDs is done by computing the maximal ID in the initial condition and then attributing two batchs of ID to each rank.
The size of each batch is computed in the same way than the extra particles in order to ensure that we will have enough ID.
Every time that a new ID is requested, the next available ID in the first batch is used.
If the last available ID in the first batch is requested, the first one is replaced by the second one and the second one is flagged as being used.
If the second batch is also fully used, the code will exit with an error message. Due to the size of the slices, the code should run out of extra particles
before reaching this point.

At each rebuild, the ranks will request a new batch if required and communicate about it.
