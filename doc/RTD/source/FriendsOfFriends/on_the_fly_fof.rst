.. Friends Of Friends
   Matthieu Schaller 15th June 2019

.. _fof_on_the_fly_label:

On-the-fly Friends-Of-Friends and Black Holes seeding
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The main purpose of the on-the-fly FOF calls is to identify haloes during a
cosmological simulation in order to seed some of them with black holes
based on physical considerations.

**In this mode, no group catalog is written to the disk. The resulting list
of haloes is only used internally by SWIFT.**

Once the haloes have been identified by the FOF code, SWIFT will iterate
over the list of groups and will check whether each halo obeys the
following criteria:

  * Is above a user-specified mass threshold (typically
    :math:`10^{10}\rm{M}_\odot` or more).
  * Contains *at least* one gas particle.
  * Does *not* contain any already existing black hole particle.

If a group satisfies all these requirements then the *densest* gas particle
(based on their SPH density) in the group will be converted to a black
hole. In practice this means that the gas particle is removed and a black
hole particle with the same dynamical mass is inserted at the same position
in the domain and with the same velocity. Additional properties are copied
from the gas particle to the black hole particle depending on the specifics
of the sub-grid model in use.
