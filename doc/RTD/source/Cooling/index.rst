.. Equation of State
   Loic Hausammann, 7th April 2018


Cooling
=======

Currently, we have 5 different cooling (EAGLE, Grackle, const-lambda, const-du and none).
Three of them are easily solved analytically (const-lambda, const-du and none) while the two last requires complex chemical networks.


Equations
---------

The first table compares the different analytical cooling while the next ones are specific to a given cooling.
The quantities are the internal energy (\\( u \\)), the density \\( rho \\), the element mass fraction (\\( X_i \\)), the cooling function (\\(\\Lambda\\), the proton mass (\\( m_H \\)) and the time step condition (\\( t\_\\text{step}\\)).
If not specified otherwise, all cooling contains a temperature floor avoiding negative temperature.

.. csv-table:: Analytical Cooling
   :header: "Variable", "Const-Lambda", "Const-du", "None"

   "\\( \\frac{ \\mathrm{d}u }{ \\mathrm{d}t } \\)", "\\( -\\Lambda \\frac{\\rho^2 X_H^2}{\\rho m_H^2} \\)", "const", "0"
   "\\( \\Delta t\_\\text{max} \\)", "\\( t\_\\text{step} \\frac{u}{\\left|\\frac{ \\mathrm{d}u }{ \\mathrm{d}t }\\right|} \\)", "\\( t\_\\text{step} \\frac{u}{\\ \\left| \\frac{ \\mathrm{d}u }{ \\mathrm{d}t }\\right|} \\)", "None"


Grackle
~~~~~~~
   
Grackle is a chemistry and cooling library presented in B. Smith et al. 2016 (do not forget to cite if used).
Four different modes are available: equilibrium, 6 species network (H, H\\( ^+ \\), e\\( ^- \\), He, He\\( ^+ \\) and He\\( ^{++} \\)), 9 species network (adds H\\(^-\\), H\\(_2\\) and H\\(_2^+\\)) and 12 species (adds D, D\\(^+\\) and HD).
It also includes some self-shielding methods and UV background.
In order to use the Grackle cooling, you will need to provide an HDF5 table computed by Cloudy.

Eagle
~~~~~

TODO

How to Implement a New One
--------------------------

We recommend to split the cooling in two functions: a first one computing the cooling rate and the second one that apply 

