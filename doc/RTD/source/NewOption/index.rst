.. Equation of State
   Loic Hausammann, 7th April 2018

.. _new_option:
   
Adding a New Option
===================

The following steps are required for any new options (such as new :ref:`hydro`, :ref:`chemistry`, :ref:`cooling`, :ref:`equation_of_state`, :ref:`stars` or :ref:`gravity`)
   
In order to add a new option
1. Create a new subdirectory inside the option directory (e.g. `src/equation_of_state` or `src/hydro`) with an explicit name.

2. Create the required new files (depending on your option, you will need different files).
   Copy the structure of the most simple option (e.g. `src/hydro/Gadget2`, `src/gravity/Default`, `src/stars/Default`, `src/cooling/none`, `src/chemistry/none` or `src/equation_of_state/ideal_gas`)

3. Add the right includes in the option file (e.g. `src/hydro.h`, `src/gravity.h`, `src/stars.h`, `src/cooling.h`, `src/chemistry.h` or `src/equation_of_state.h`) and the corresponding io file if present.

4. Add the new option in `configure.ac`.
   This file generates the `configure` script and you just need to add a new option under the right `case`.

5. Add your files in `src/Makefile.am`.
   In order to generate the Makefiles during the configuration step, a list of file is required.
   In `nobase_noinst_HEADERS`, add your new header files.

6. Update the documentation.
   Add your equations/documentation in `doc/RTD`.
