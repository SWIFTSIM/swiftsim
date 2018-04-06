.. Equation of State
   Loic Hausammann, 6th April 2018


Equation of State
=================

Currently (if the documentation was well updated), we have two different equation of states implemented: ideal gas and isothermal.
They describe the relations between our main thermodynamical variables: the internal energy (\\(u\\)), the density (\\(\\rho\\)), the entropy (\\(A\\)) and the pressure (\\(P\\)).

Equations
---------
In the following section, the variables not yet defined are: \\(\\gamma\\) for the adiabatic index and \\( c_s \\) for the speed of sound.

.. csv-table:: Ideal Gas
   :header: "Variable", "A", "u", "P"
	   
   "A", "", "\\( \\left( \\gamma - 1 \\right) u \\rho^{1-\\gamma} \\)", "\\(P \\rho^{-\\gamma} \\)"
   "u", "\\( A \\frac{ \\rho^{ \\gamma - 1 } }{\\gamma - 1 } \\)", "", "\\(\\frac{1}{\\gamma - 1} \\frac{P}{\\rho}\\)"
   "P", "\\( A \\rho^\\gamma \\)", "\\( \\left( \\gamma - 1\\right) u \\rho \\)", ""
   "\\(c_s\\)", "\\(\\sqrt{ \\gamma \\rho^{\\gamma - 1} A}\\)", "\\(\\sqrt{ u \\gamma \\left( \\gamma - 1 \\right) } \\)", "\\(\\sqrt{ \\frac{\\gamma P}{\\rho} }\\)"


.. csv-table:: Isothermal Gas
   :header: "Variable", "A", "u", "P"

	    
   "A", "", "\\(\\left( \\gamma - 1 \\right) u \\rho^{1-\\gamma}\\)", "" 
   "u", "", "const", ""
   "P", "", "\\(\\left( \\gamma - 1\\right) u \\rho \\)", ""
   "\\( c_s\\)", "", "\\(\\sqrt{ u \\gamma \\left( \\gamma - 1 \\right) } \\)", ""


How to Implement a New One
--------------------------

If you wish to implement a new equation of state, you will need to do a few different things:

1. Create a new subdirectory inside `src/equation_of_state` (with the same name than your new option).

2. Create a new file containing your equation of state.
   It should contains the definition of `eos_parameters`, IO functions and transformations between the different variables: \\(u(\\rho, A)\\), \\(u(\\rho, P)\\), \\(P(\\rho,A)\\), \\(P(\\rho, u)\\), \\(A(\\rho, P)\\), \\(A(\\rho, u)\\), \\(c_s(\\rho, A)\\), \\(c_s(\\rho, u)\\) and \\(c_s(\\rho, P)\\). See other equation of state files to have implementation details.

3. Add the right include in `src/equation_of_state.h`.
   In this file, we are selecting the equation of state from the configuration file.

4. Add the new option in `configure.ac`.
   This file generates the `configure` script.
   Simply add the right option name (below `case "$with_eos" in`) and definition (using `AC_DEFINE`) in the equation of state switch.

5. Add your files in `src/Makefile.am`.
   In order to generate the Makefiles during the configuration step, a list of file is required.
   In `nobase_noinst_HEADERS`, add your new files.

6. Update the documentation.
   Add your equations in `doc/RTD/source/EquationOfState/index.rst`.
