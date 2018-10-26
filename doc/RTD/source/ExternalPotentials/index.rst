.. External potentials in SWIFT
   Folkert Nobels, 25th October 2018
   
External Potentials in SWIFT
============================

SWIFT can be run with an external potential on this page we will summaraize the
current potentials which can be run with SWIFT and how to implement your own 
potential in SWIFT.

Current External Potentials
---------------------------


How to implement your own potential
-----------------------------------

The first step in implementing your own potential is making a directory of your
potential in the ``src/potential`` folder and creating a file in the folder 
called ``potential.h``.

Configuring the potential 
^^^^^^^^^^^^^^^^^^^^^^^^^

THe first thing is copying a ``potential.h`` file from an already implemented 
potential. In this potential the ``ifdef`` statements need to be changed to the
specific potential and the ``struct`` and ``potential_init_backend`` need to be
changed such that it uses your potential and reads the correct potential from
the parameter file during running the program.

Add the potential to the ``potential.h`` file in the ``src`` directory.

Add the potential to the ``config.h`` header such that we can configure the 
program to use the correct potential.
