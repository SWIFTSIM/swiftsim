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

CHIMES also requires the CVODE library of Ordinary Differential Equation (ODE) solvers, from the Sundials package (v4.0 or later). If Sundials is not already installed on your system, it can be downloaded `here <https://computing.llnl.gov/projects/sundials/sundials-software>`_. 

Once you have downloaded Sundials from the above website, you can then untar it and build it. You will need ``cmake`` to build Sundials. The INSTALL_GUIDE.pdf included with the Sundials download describes this process in detail, but in brief these are the steps that you will need: 

.. code-block:: bash

  tar -zxvf sundials-5.1.0.tar.gz 
  cd sundials-5.1.0 
  mkdir build 
  cd build 
  cmake -DCMAKE_INSTALL_PREFIX=/path/to/install/dir/ -DBUILD_ARKODE=OFF -DBUILD_CVODE=ON -DBUILD_CVODES=OFF -DBUILD_IDA=OFF -DBUILD_IDAS=OFF -DBUILD_KINSOL=OFF -DBUILD_SHARED_LIBS=ON -DBUILD_STATIC_LIBS=ON -DCMAKE_C_FLAGS="-O2" -DEXAMPLES_ENABLE_C=OFF ../
  make
  make install

In the above example, after untar'ing the sundials download we create a new ``build`` directory inside ``sundials-5.1.0``, and then ``cd build`` before running ``cmake``. This is needed because the build directory has to be different from the source directory. 

Also, we only need the CVODE library, so I have switched off building the various other libraries that are included in the Sundials package, as these are not required by CHIMES. 

Finally, you will need to set the ``-DCMAKE_INSTALL_PREFIX`` path to a directory where you have write access. This is where the libraries will be installed. 

Once you have installed Sundials, make a note of the directory that it has been installed in, as you may need to pass this to the configure script when you build SWIFT (see the next section). You will also need to add ``/path/to/install/dir/lib`` (or possibly ``/path/to/install/dir/lib64``, depending on your system) to your `LD_LIBRARY_PATH` environment variable, both when you compile SWIFT and when you run it. 

