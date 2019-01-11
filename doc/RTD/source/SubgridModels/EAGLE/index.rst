.. EAGLE sub-grid model
   Matthieu Schaller, 20th December 2018


EAGLE model
===========

This section of the documentation gives a brief description of the
different components of the EAGLE sub-grid model. We mostly focus on
the parameters and values output in the snapshots.

Chemical tracers
~~~~~~~~~~~~~~~~

The gas particles in the EAGLE model carry metal abundance information
in the form of metal mass fractions. We follow the following 9
elements: `H`, `He`, `C`, `N`, `O`, `Ne`, `Mg`, `Si` and `Fe`. We
additionally follow the total metal mass fraction (i.e. absolute
metallicity) `Z`. This is typically larger than the sum of the 7
metals that are individually traced since this will also contain the
contribution of all the elements that are not individually followed.
We note that all of definitions are independent of any definition of
solar the solar metallicity :math:`Z_\odot` or of any solar abundance
pattern.

As part of the diagnostics, we additionally trace the elements coming
from the different stellar evolution channels. We store for each
particle the total mass coming from all the SNIa that enriched that
particle and the metal mass fraction from SNIa. This is the fraction
of the *total* gas mass that is in the form of metals originating from
SNIa stars. By construction this fraction will be smaller than the
total metal mass fraction. The same tracers exist for the SNII and AGB
channels. Finally, we also compute the iron gas fraction from
SNIa. This it the fraction of the *total* gas mass that is made of
iron originating from SNIa explosions. 

We finally also compute the smoothed versions of the individual
element mass fractions, of the total metal mass fractions, and of the
iron gas fraction from SNIa.

The chemistry module in ``src/chemistry/EAGLE`` includes all the arrays
that are added to the particles and the functions used to compute the
smoothed elements.

When a star is formed (see below), it inherits all the chemical
tracers of its parent gas particle.

In the snapshots, we output for each gas and star particle:

+----------------------------------+-------------------------------------+-----------+-----------------------------+
| Name                             | Description                         | Units     | Comments                    |
+==================================+=====================================+===========+=============================+
| ``ElementAbundance``             | | Fraction of the gas/star mass     | [-]       | | Array of length           |
|                                  | | in the different elements         |           | | 9 for each particle       |
+----------------------------------+-------------------------------------+-----------+-----------------------------+
| ``SmoothedElementAbundance``     | | Fraction of the gas/star mass     | [-]       | | Array of length           |
|                                  | | in the different elements         |           | | 9 for each particle       |
|                                  | | smoothed over SPH neighbours      |           |                             |
+----------------------------------+-------------------------------------+-----------+-----------------------------+
| ``Metallicity``                  | | Fraction of the gas/star mass     | [-]       |                             |
|                                  | | in *all* metals                   |           |                             |
+----------------------------------+-------------------------------------+-----------+-----------------------------+
| ``SmoothedMetallicity``          | | Fraction of the gas/star mass     | [-]       |                             |
|                                  | | in *all* metals                   |           |                             |
|                                  | | smoothed over SPH neighbours      |           |                             |
+----------------------------------+-------------------------------------+-----------+-----------------------------+
| ``TotalMassFromSNIa``            | | Total mass of the gas/star        | [U_M]     |                             |
|                                  | | that was produced by enrichment   |           |                             |
|                                  | | from SNIa stars                   |           |                             |
+----------------------------------+-------------------------------------+-----------+-----------------------------+
| ``MetalMassFracFromSNIa``        | | Fraction of the *total* gas/star  | [-]       |                             |
|                                  | | mass that is in metals produced   |           |                             |
|                                  | | by enrichment from SNIa stars     |           |                             |
+----------------------------------+-------------------------------------+-----------+-----------------------------+
| ``TotalMassFromAGB``             | | Total mass of the gas/star        | [U_M]     |                             |
|                                  | | that was produced by enrichment   |           |                             |
|                                  | | from AGB stars                    |           |                             |
+----------------------------------+-------------------------------------+-----------+-----------------------------+
| ``MetalMassFracFromAGB``         | | Fraction of the *total* gas/star  | [-]       |                             |
|                                  | | mass that is in metals produced   |           |                             |
|                                  | | by enrichment from AGB star       |           |                             |
+----------------------------------+-------------------------------------+-----------+-----------------------------+
| ``TotalMassFromSNII``            | | Total mass of the gas/star        | [U_M]     |                             |
|                                  | | that was produced by enrichment   |           |                             |
|                                  | | from SNII stars                   |           |                             |
+----------------------------------+-------------------------------------+-----------+-----------------------------+
| ``MetalMassFracFromSNII``        | | Fraction of the gas/star mass     | [-]       |                             |
|                                  | | that is in metals produced by     |           |                             |
|                                  | | enrichment from SNII stars        |           |                             |
+----------------------------------+-------------------------------------+-----------+-----------------------------+
| ``IronMassFracFromSNIa``         | | Fraction of the *total* gas/star  | [-]       |                             |
|                                  | | mass in *iron* produced produced  |           |                             |
|                                  | | by enrichment from SNIa stars     |           |                             |
+----------------------------------+-------------------------------------+-----------+-----------------------------+
| ``SmoothedIronMassFracFromSNIa`` | | Fraction of the *total* gas/star  | [-]       |                             |
|                                  | | mass in *iron* produced produced  |           |                             |
|                                  | | by enrichment from SNIa stars     |           |                             |
|                                  | | smoothed over SPH neighbours      |           |                             |
+----------------------------------+-------------------------------------+-----------+-----------------------------+

The stars will lose mass over their lifetime (up to ~45%). The
fractions will remain unchanged but if one is interested in computing
an absolute metal mass (say) for a star, the ``InitialMass`` (see
below) of the star must be used.

The chemistry model only requires a small number of parameters to be
specified in the `EAGLEChemistry` section of the YAML file. These are
the initial values of the metallicity and element mass
fractions. These are then applied at the start of a simulation to
*all* the *gas* particles. All 9 elements have to be specified An
example section, for primordial abundances (typical for a cosmological
run), is:

.. code:: YAML

   EAGLEChemistry:
     Metallicity:                0.    # Mass fraction in all metals
     InitAbundance_Hydrogen:     0.755 # Mass fraction in Hydrogen
     InitAbundance_Helium:       0.245 # Mass fraction in Helium
     InitAbundance_Carbon:       0.    # Mass fraction in Carbon
     InitAbundance_Nitrogen:     0.    # Mass fraction in Nitrogen
     InitAbundance_Oxygen:       0.    # Mass fraction in Oxygen
     InitAbundance_Neon:         0.    # Mass fraction in Neon
     InitAbundance_Magnesium:    0.    # Mass fraction in Magnesium
     InitAbundance_Silicon:      0.    # Mass fraction in Silicon
     InitAbundance_Iron:         0.    # Mass fraction in Iron

Whilst one would use the following values for solar abundances
(typical for an idealised low-redshift run):

.. code:: YAML

   EAGLEChemistry:
     Metallicity:                0.014        # Mass fraction in all metals
     InitAbundance_Hydrogen:     0.70649785   # Mass fraction in Hydrogen
     InitAbundance_Helium:       0.28055534   # Mass fraction in Helium
     InitAbundance_Carbon:       2.0665436e-3 # Mass fraction in Carbon
     InitAbundance_Nitrogen:     8.3562563e-4 # Mass fraction in Nitrogen
     InitAbundance_Oxygen:       5.4926244e-3 # Mass fraction in Oxygen
     InitAbundance_Neon:         1.4144605e-3 # Mass fraction in Neon
     InitAbundance_Magnesium:    5.907064e-4  # Mass fraction in Magnesium
     InitAbundance_Silicon:      6.825874e-4  # Mass fraction in Silicon
     InitAbundance_Iron:         1.1032152e-3 # Mass fraction in Iron


     
Gas cooling: Wiersma+2009a
~~~~~~~~~~~~~~~~~~~~~~~~~~

The gas cooling is based on the redshift-dependent tables of `Wiersma et
al. (2009) <http://adsabs.harvard.edu/abs/2009MNRAS.393...99W>`_ that include
element-by-element cooling rates for the 11 elements (`H`, `He`, `C`, `N`, `O`,
`Ne`, `Mg`, `Si`, `S`, `Ca` and `Fe`) that dominate the total rates. The tables
assume that the gas is in ionization equilibrium with the cosmic microwave
background (CMB) as well as with the evolving X-ray and UV background from
galaxies and quasars described by the model of `Haardt & Madau (2001)
<http://adsabs.harvard.edu/abs/2001cghr.confE..64H>`_. Note that this model
ignores *local* sources of ionization, self-shielding and non-equilibrium
cooling/heating. The tables can be obtained from this `link
<http://virgodb.cosma.dur.ac.uk/swift-webstorage/CoolingTables/EAGLE/coolingtables.tar.gz>`_
which is a re-packaged version of the `original tables
<http://www.strw.leidenuniv.nl/WSS08/>`_

The Wiersma tables containing the cooling rates as a function of redshift,
Hydrogen number density, Helium fraction (:math:`X_{He} / (X_{He} + X_{H})`) and
element abundance relative to the solar abundance pattern assumed by the tables
(see equation 4 in the original paper). As the particles do not carry the mass
fraction of `S` and `Ca`, we compute the contribution to the cooling rate of
these elements from the abundance of `Si`. More specifically, we assume that
their abundance by mass relative to the table's solar abundance pattern is the
same as the relative abundance of `Si` (i.e. :math:`[Ca/Si] = 0` and
:math:`[S/Si] = 0`). Users can optionally modify the ratios used for `S` and
`Ca`.

Above the redshift of Hydrogen re-ionization we use the extra table containing
net cooling rates for gas exposed to the CMB and a UV + X-ray background at
redshift nine truncated above 1 Rydberg. At the redshift or re-ionization, we
additionally inject a fixed user-defined amount of energy per unit mass to all
the gas particles.

In addition to the tables we inject extra energy from Helium re-ionization using
a Gaussian model with a user-defined redshift for the centre, width and total
amount of energy injected per unit mass.

For non-cosmological run, we use the :math:`z = 0` table and the interpolation
along the redshift dimension then becomes a trivial operation.

The cooling itself is performed using an implicit scheme (see the theory
documents) which for small values of the cooling rates is solved explicitly. For
larger values we use a bisection scheme. Users can alternatively use a
Newton-Raphson method that in some cases runs faster than the bisection
method. If the Newton-Raphson method does not converge after a few steps, the
code reverts to a bisection scheme, that is guaranteed to converge. The cooling
rate is added to the calculated change in energy over time from the other
dynamical equations. This is different from other commonly used codes in the
literature where the cooling is done instantaneously.

We note that the EAGLE cooling model does not impose any restriction on the
particles' individual time-steps. The cooling takes place over the time span
given by the other conditions (e.g the Courant condition).

The cooling model is driven by a small number of parameter files in the
`EAGLECooling` section of the YAML file. These are the re-ionization parameters,
the path to the tables and optionally the modified abundances of `Ca` and `S` as
well as the flag to attempt using the Newton-Raphson scheme to solve the
implicit problem. A valid section of the YAML file looks like:

.. code:: YAML

   EAGLECooling:
     dir_name:     /path/to/the/Wiersma/tables/directory # Absolute or relative path
     H_reion_z:            11.5      # Redhift of Hydrogen re-ionization
     He_reion_z_centre:     3.5      # Centre of the Gaussian used for Helium re-ionization
     He_reion_z_sigma:      0.5      # Width of the Gaussian used for Helium re-ionization
     He_reion_eV_p_H:       2.0      # Energy injected in eV per Hydrogen atom for Helium re-ionization.

And the optional parameters are:

.. code:: YAML

   EAGLECooling:
     Ca_over_Si_in_solar:       1.0 # (Optional) Value of the Calcium mass abundance ratio to solar in units of the Silicon ratio to solar. Default value: 1.
     S_over_Si_in_solar:        1.0 # (Optional) Value of the Sulphur mass abundance ratio to solar in units of the Silicon ratio to solar. Default value: 1.
     newton_integration:        0   # (Optional) Set to 1 to use the Newton-Raphson scheme for the explicit cooling problem.



Particle tracers
~~~~~~~~~~~~~~~~

Star formation: Schaye+2008
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Stellar enrichment: Wiersma+2009b
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Supernova feedback: Dalla Vecchia+2012
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Black-hole creation
~~~~~~~~~~~~~~~~~~~

Black-hole accretion
~~~~~~~~~~~~~~~~~~~~

AGN feedback
~~~~~~~~~~~~
