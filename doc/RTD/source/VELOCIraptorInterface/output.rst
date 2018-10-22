.. VELOCIraptor output
   Folkert Nobels 12th of October

VELOCIraptor Output
===================

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: Contents: 

In general VELOCIraptor outputs six files per snapshot, of which 2 files are
for unbound particles specifically.  In this part we will explain what is
inside the different files.

Catalog_groups file
-------------------

The first file that is output by VELOCIraptor is the ``.catalog_group`` file,
this file contains all the information that is group specific, the interesting
data in the ``.catalog_group`` files are: 

+ The ``group_size``: gives a list of all the halos and the number of particles
  in the halo, this list is numbered from 0 until the number of groups minus
  one. It is important that the groups are not ordered in any way [#order]_ 
+ The ``Num_of_groups`` or ``Total_num_of_groups``: gives the total number of
  groups in the snapshot.
+ The ``Offset`` list: This list gives the offset off the particles. In the
  output of VELOCIraptor there is no file which has an ID for every particle
  and a corresponding group, rather the particles are ordered according to in
  which group they are. So if we want to access the particles in group 0, we
  need to look at the particles from ``Offset[0]`` until ``Offset[1]`` in the
  ``.catalog_particles`` hdf5 file. In general this means that for group N we
  need to look at particles ``Offset[N]`` until ``Offset[N+1]``. 
+ The ``Offset_unbound`` list: This list works exactly the same as the
  ``Offset`` list only this list is for the gravitational unbound particles.

Catalog_particles file
----------------------

The second file that is produced by VELOCIraptor is the ``.catalog_particles``
file, this file contains mainly all the IDs of the particles and mainly has two
interesting things:

+ The ``Num_of_particles_in_groups`` and ``Num_of_particles_in_groups``
  parameter: Gives the total number of particles in the file or which are found
  in the halo. 
+ The ``Particle_IDs``: The list of particles as sorted by halo, in which halo
  the individual particles are present can be found by using the
  ``.catalog_group`` file and the corresponding ``Offset`` list. 

Besides the ``.catalog_particles`` file, there is also a
``.catalog_particles.unbound`` file, this file contains the same information
but only for the unbound particles, a particle can only be present in one of
these two lists. 

Catalog_parttypes file
----------------------

The third file that is produced by VELOCIraptor is the ``.catalog_parttypes``
file, this file contains the information what type of particle every particle
is, ordered the same as in ``Particle_IDs`` in ``.catalog_particles``. There
are only two interesting parameters of the file which are:

+ The ``Num_of_particles_in_groups`` parameter: Gives the total number of
  particles in the file which are in a halo.
+ The ``Particle_types`` list: Gives a list of particles types similar to the
  snap shots (0 - gas, 1 - dm, 4 - stars).

Besides the ``.catalog_parttypes`` file, there is also a
``.catalog_parttypes.unbound`` file, this file contains this information for
the unbound particles.

Properties file
---------------

The Fourth file is the ``.properties`` file, this file contains mainy physical
useful information of the corresponding halos. This can be divided in several
useful groups of physical parameters, on this page we have divided the several
variables which are present in the ``.properties`` file.

Critical Density related:
"""""""""""""""""""""""""

+ ``Mass_200crit``: The mass of a halo with an over density on average of
  :math:`\Delta=200` based on the critical density of the Universe 
  (:math:`M_{200}`).
+ ``R_200crit``: The :math:`R_{200}` radius of the halo based on the 
  critical density of the Universe

Mean Density related:
"""""""""""""""""""""

+ ``Mass_200mean``: The mass of a halo with an over density on average of
  :math:`\Delta=200` based on the mean density of the Universe 
  (:math:`M_{200}`).
+ ``R_200mean``: The :math:`R_{200}` radius of the halo based on the 
  mean density ofthe Universe.

Virial properties:
""""""""""""""""""

+ ``Mvir``: The virial mass of the halos.
+ ``Rvir``: The virial radius of the halo (:math:`R_{vir}`).

NFW profile properties:
"""""""""""""""""""""""

+ ``cNFW``: The concentration of the halo.
+ ``Xc``, ``Yc`` and ``Zc``: The x,y and z centre positions of the halos 
  [#center]_.
+ ``Xc_gas``, ``Yc_gas``, ``Zc_gas``: The offset of the centre positions of
  the halo based on the gas, to find the position of the gas the offsets 
  need to be added to ``Xc``, ``Yc`` and ``Zc``. 

Velocity Dispersion related:
""""""""""""""""""""""""""""

+ The complete velocity dispersion tensor (:math:`\sigma_{ij}`) which has 
  an array per component which gives the value for all the halos. In 
  general these components are called ``veldisp_ij`` in which i and j are 
  given by ``x``, ``y`` or ``z``. This means that there are nine 
  components stored in the ``.properties`` file. This omits the fact 
  that the dispersion tensor by nature is a symmetric tensor. All the 
  components are given by: 
  ``veldisp_xx``, ``veldisp_xy``, ``veldisp_xz``, ``veldisp_yx``, 
  ``veldisp_yy``, ``veldisp_yz``, ``veldisp_zx``, ``veldisp_zy``, 
  and ``veldisp_zz`` [#velodisp]_.
+ ``sigV``, the scalar velocity dispersion which corresponds with the 
  trace of the velocity dispersion tensor 
  (:math:`\sigma = \text{Tr}(\sigma_{ij})`).


Intertia Tensor properties:
"""""""""""""""""""""""""""

+ ``eig_ij``: Are the normalized eigenvectors of the velocity 
  dispersion, :math:`\sigma_{ij}`.
+ The eigenvalue ratios: 
  1. ``q`` is the semi-major over major; 
  2. ``s`` is the minor over major.
+ The eigenvalue ratios for only the gas, similar to all particles:
  1. ``q_gas`` is the semi-major over major for only gas; 
  2. ``s_gas`` is the minor over major for only gas.

Spin parameters:
""""""""""""""""

+ ``lambda_b`` is the bullock spin parameter, see the paper by Bullock et al. (2001) [#Bullock]_. 

And some more properties which do not have a specified category yet:

  + ``Mass_FOF``: The friends-of-friends mass of the halos.

  + ``Other parameters``: Soon



.. [#order] In most cases more massive groups appear earlier in the list, but 
   this is not guaranteed for larger simulations. The order of the groups is 
   more a matter of the way that VELOCIraptor searches instead of a physical 
   reason.
.. [#center] This is not the average positions of the halos particles, but
   the halo position found by the VELOCIraptor algorithm. This includes a 
   fit for all the parameters including the gas particles or other types of
   particles.
.. [#velodisp] In the velocity dispersion tensor ( :math:`\sigma_{ij}` )  
   the following relations are satisfied between components:

   + :math:`\sigma_{xy}=\sigma_{yx}`
   + :math:`\sigma_{xz}=\sigma_{zx}`
   + :math:`\sigma_{yz}=\sigma_{yz}`
.. [#Bullock] The Bullock spin parameter is given by 
   :math:`\lambda = \frac{J}{\sqrt{2}MVR}`, for more information see 
   https://arxiv.org/abs/astro-ph/0011001. 
