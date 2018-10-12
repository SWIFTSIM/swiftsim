.. VELOCIraptor stand alone 
   Folkert Nobels 12th October 2018

Stand alone VELOCIraptor configuration
======================================


.. toctree::    
   :maxdepth: 2    
   :hidden:    
   :caption: Contents: 
   
Besides running VELOCIraptor on the fly when using SWIFT, it is also possible
to run VELOCIraptor alone without using SWIFT. In this section we explain how 
VELOCIraptor can be run stand alone without using SWIFT.

Setting up VELOCIraptor
-----------------------

The first step is setting up VELOCIraptor, this requires us to download the 
git repository as::
  
  git clone https://github.com/pelahi/VELOCIraptor-STF

Similar to the SWIFT with VELOCIraptor configuration, we can use the 
swift-interface branch to analyze individual snapshots. We can use this branch
by doing::

  cd VELOCIraptor-STF
  git fetch
  git checkout swift-interface

Again we need to copy the default SWIFT config file to a other config file by
doing::

  cd stf
  cp Makefile.config.SWIFT-template Makefile.config

Similar to configuring VELOCIraptor with swift we need to change the first 20
lines of ``Makefile.config`` to work with our compiler, but we also need to 
change the fact that we do not use the swift-interface but the standalone 
version of the code, so change ``SWIFTINTERFACE="on"`` to 
``SWIFTINTERFACE="off"``.

Compiling VELOCIraptor
----------------------

Compoling goes completely different as compared to the on the fly halo finder
configuration with SWIFT. In this case we can compile the code as::

  make 

After this an additional folder is created in ``VELOCIraptor-stf/stf`` called
``bin``, in which the binary files of ``stf-gas`` is present (assuming you 
run a simulation with SPH [#nosph]_)

Running VELOCIraptor on a Snapshot
----------------------------------

After the code is compile the next step is calculating 



.. [#nosph] In the case that in the ``Makefile.config`` it is indicate that the 
   simulation does only contain dark matter this will reflect back on the 
   generated binary file. So ``stf-gas`` will change to ``stf`` in the case of 
   a dark matter only simulation.

