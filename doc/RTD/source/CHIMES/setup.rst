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

Once you have downloaded Sundials from the above website, you can then untar it and build it. You will need ``cmake`` to build Sundials. The INSTALL_GUIDE.pdf included with the Sundials download describes this process in detail, but in brief these are the steps that you will need: 

.. code-block:: bash

  tar -zxvf sundials-2.6.0.tar.gz 
  cd sundials-2.6.0 
  mkdir build 
  cd build 
  cmake -DCMAKE_INSTALL_PREFIX=/path/to/install/dir/ -DBUILD_ARKODE=OFF -DBUILD_CVODE=ON -DBUILD_CVODES=OFF -DBUILD_IDA=OFF -DBUILD_IDAS=OFF -DBUILD_KINSOL=OFF -DBUILD_SHARED_LIBS=ON -DBUILD_STATIC_LIBS=OFF -DCMAKE_C_FLAGS="-O2" -DEXAMPLES_ENABLE=OFF ../
  make
  make install

In the above example, after untar'ing the sundials download we create a new ``build`` directory inside ``sundials-2.6.0``, and then ``cd build`` before running ``cmake``. This is needed because the build directory has to be different from the source directory. 

Also, we only need the CVODE library, so I have switched off building the various other libraries that are included in the Sundials package, as these are not required by CHIMES. 

Finally, you will need to set the ``-DCMAKE_INSTALL_PREFIX`` path to a directory where you have write access. This is where the libraries will be installed. 

Once you have installed Sundials, make a note of the directory that it has been installed in, as you will need to pass this to the configure script when you build SWIFT (see the next section). 

