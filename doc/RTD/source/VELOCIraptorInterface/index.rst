.. VELOCIraptor Interface
   Folkert Nobels, 8th October 2018

VELOCIraptor Interface
======================

In SWIFT it is possible to run a cosmological simulation and at the same time do an on the fly halo finder at specific predefined intervals. 
Because of this we will explain on this page how we can set up an simulation using VELOCIraptor (formerly STructure Finder). 
After this we will explain what the outputs of VELOCIraptor will be.

Configuring SWIFT
-----------------

Setting up VELOCIraptor
~~~~~~~~~~~~~~~~~~~~~~~

Before we can run SWIFT with VELOCIraptor we first need to download VELOCIraptor. This can be done by cloning the repository on GitHub_::
  git clone https://github.com/pelahi/VELOCIraptor-STF

Currently the best version that works with SWIFT is the swift-interface branch of VELOCIraptor, to get this branch use::
  cd VELOCIraptor-STF
  git fetch
  git checkout swift-interface

To get the default that works with SWIFT simply copy the SWIFT template file in the ``Makefile.config``::
  cp Makefile.config.SWIFT-template Makefile.config

Depending on your compiler you want to change the first 20 lines of your ``Makefile.config`` to work with your compiler and either have MPI on or off. 


Compiling VELOCIraptor
~~~~~~~~~~~~~~~~~~~~~~

The next part will be to compile VELOCIraptor, this can simply be done by using::
  make 

After the compilation of your code, there is an additional folder created in the ``VELOCIraptor-stf/stf`` directory called ``lib`` this directory has the libary of VELOCIraptor and is required to run SWIFT with VELOCIraptor.

Compiling SWIFT
~~~~~~~~~~~~~~~
The next part is compiling SWIFT with VELOCIraptor and assumes you already downloaded SWIFT from the GitLab_, this can be done by running::
  ./autogen.sh
  ./configure --with-velociraptor=/path/to/VELOCIraptor-STF/stf/lib
  make 

In which ``./autogen.sh`` only needs to be run once after the code is cloned from the GitLab_, and ``/path/to/`` is the path to the ``VELOCIraptor-STF`` directory on your machine. In general ``./configure`` can be run with other options as desired. After this we can run SWIFT with VELOCIraptor, in the case of the Small Cosmological Volume example this will can be run as::
  cd examples/SmallCosmoVolume 
  ../swift -c -s -G -x -t 8 small_cosmo_volume.yml
In which ther is an additional ``-x`` option which activates the VELOCIraptor interface. 


VELOCIraptor Output
-------------------

In general VELOCIraptor outputs six files per snapshot, of which 2 files are for unbound particles specifically. 
In this part we will explain what is inside the different files.

Catalog_groups file
~~~~~~~~~~~~~~~~~~~

The first file that is output by VELOCIraptor is the ``.catalog_group`` file, this file contains all the information that is group specific, the interesting data in the ``.catalog_group`` files are: 
- The ``group_size``: gives a list of all the halos and the number of particles in the halo, this list is numbered from 0 until the number of groups minus one. 
- The ``Num_of_groups`` or ``Total_num_of_groups``: gives the total number of groups in the snapshot.
- The ``Offset`` list: This list gives the offset off the particles. In the output of VELOCIraptor there is no file which has an ID for every particle and a corresponding group, rather the particles are ordered according to in which group they are. So if we want to access the particles in group 0, we need to look at the particles from ``Offset[0]`` until ``Offset[1]`` in the ``.catalog_particles`` hdf5 file. In general this means that for group N we need to look at particles ``Offset[N]`` until ``Offset[N+1]``. 
- The ``Offset_unbound`` list: This list works exactly the same as the ``Offset`` list only this list is for the gravitational unbound particles.

Catalog_particles file
~~~~~~~~~~~~~~~~~~~~~~

The second file that is produced by VELOCIraptor is the ``.catalog_particles`` file, this file contains mainly all the IDs of the particles and mainly has two interesting things:
- The ``Num_of_particles_in_groups`` and ``Num_of_particles_in_groups`` parameter: Gives the total number of particles in the file or which are found in the halo. 
- The ``Particle_IDs``: The list of particles as sorted by halo, in which halo the individual particles are present can be found by using the ``.catalog_group`` file and the corresponding ``Offset`` list. 

Besides the ``.catalog_particles`` file, there is also a ``.catalog_particles.unbound`` file, this file contains the same information but only for the unbound particles, a particle can only be present in one of these two lists. <!-- check -->

Catalog_parttypes file
~~~~~~~~~~~~~~~~~~~~~~

The third file that is produced by VELOCIraptor is the ``.catalog_parttypes`` file, this file contains the information what type of particle every particle is, ordered the same as in ``Particle_IDs`` in ``.catalog_particles``. There are only two interesting parameters of the file which are:
- The ``Num_of_particles_in_groups`` parameter: Gives the total number of particles in the file which are in a halo.
- The ``Particle_types`` list: Gives a list of particles types similar to the snap shots (0 - gas, 1 - dm, 4 - stars).

Besides the ``.catalog_parttypes`` file, there is also a ``.catalog_parttypes.unbound`` file, this file contains this information for the unbound particles.

Properties file
~~~~~~~~~~~~~~~
The Fourth file is the ``.properties`` file, this file contains mainly physical usefull information of the corresponding halos. Some usefull physical parameters are:
- 




.. _GitHub: https://github.com/pelahi/VELOCIraptor-STF
.. _GitLab: https://gitlab.cosma.dur.ac.uk/swift/swiftsim
   
