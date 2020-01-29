.. CHIMES setup 
   Alexander Richings 28th January 2020 

.. _CHIMES_setup:

CHIMES setup
------------

Before configuring SWIFT with the CHIMES module enabled, you will need to download the CHIMES data files, and install the CVODE library from the Sundials package. These are described further below. 


CHIMES data files
^^^^^^^^^^^^^^^^^

To run CHIMES, you will need the CHIMES data files that contain the various reaction rate coefficients, photoionisation cross sections etc. These can be downloaded from the `chimes-data <https://bitbucket.org/richings/chimes-data>`_ Bitbucket repository. 

Make sure to use the master branch of the chimes-data repository. The legacy branch was for an old version of CHIMES, and is incompatible with the version of CHIMES that has been implemented in SWIFT. 

You will need to make a note of the path to your local chimes-data repository, as you will need to pass this as a parameter to SWIFT (see the CHIMES parameters section for more details). 


Sundials library
^^^^^^^^^^^^^^^^

CHIMES also requires the CVODE library of Ordinary Differential Equation (ODE) solvers, from the Sundials package. Sundials can be downloaded `here <https://computing.llnl.gov/projects/sundials/sundials-software>`_. Note that CHIMES has not been tested with the most recent version of CVODE. You should only use version 2.6.0 of Sundials (on the above page, scroll down to the Archive section, where it lists the old versions of Sundials). 

Once you have downloaded Sundials from the above website, you can then untar it and build it. You will need to build it as a shared library, i.e. using the -fpic option. If you don't have super user access on the machine that you are installing Sundials on, you can specify a directory where you do have write access to install it, using the --prefix=/path/to/dir option in the configure script. 

Also, note that choosing an appropriate compiler and optimisation flags when building Sundials can affect the speed of CHIMES by up to a factor ~2. If you have access to the Intel compilers, I find the following works best (but the best setup will depend on the machine that you are running on): 

.. code-block:: bash

  ./configure --prefix=/path/to/dir CC=icc CFLAGS="-O2 -fpic" F77=ifort FFLAGS="-O2 -fpic"
  make
  make install

Once you have installed Sundials, make a note of the directory that it has been installed in, as you will need to pass this to the configure script when you build SWIFT (see the next section). 
