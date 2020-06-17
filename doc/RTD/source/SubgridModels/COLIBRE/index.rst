.. COLIBRE sub-grid model
   Matthieu Schaller, 20th December 2018
   Folkert Nobels, 3th of June 2019


COLIBRE model
=============

This section of the documentation gives a brief description of the
different components of the COLIBRE sub-grid model. 

.. _COLIBRE_star_formation:

Star formation: Schmidt law
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The star formation model in the COLIBRE subgrid model consists of 3 components,
the first component are the star formation criteria, this determines which gas 
is allowed to form stars. The second component is the Schmidt law which sets 
the star formation rate of the gas according to a Schmidt law. The last 
component is the pressure law which sets the star formation rate according to
the pressure of the gas. Component 2 and 3 can never be used at the same time.
In the following part we will explain how to use the different parts of the 
star formation model.

Star forming criteria
^^^^^^^^^^^^^^^^^^^^^

In our simulations not all the gas is star forming, because of that we need to
determine in the simulation which gas is star forming. In the COLIBRE star 
formation model there are several star formation criteria that can be used to 
determine if gas is star forming. 

The first star formation criterion is an 
overdensity criterion. In the code the overdensity criterion is called 
:code:`min_over_density`, this is a mandatory criterion and we often set this
to 57.7. This criterion always needs to be satisfied in order of gas to be 
star forming (in non-cosmological runs it is satisfied for all gas). Besides 
this there are 4 additional criteria that are on their own plus the overdensity
criterion, if one of these criteria is satisfied and the overdensity criterion is 
satisfied, this is enough to allow the gas to be star forming. In cosmologial
runs the overdensity criteria is important because it prevents gas to be star
forming at high redshifts in non-collapsed objects.

The second criterion is the virial criterion, which uses the 3D velocity 
dispersion, the kernel mass and the density to determine if gas is 
gravitationally bound. In the parameter file the variable is called
:code:`alpha_virial`, this variable is mandatory and often set to 1.0. If
it is desired to not use the virial criteria a value of 0 will turn it off.

The third criterion is the subgrid temperature threshold, which allows gas to
be star forming when the subgrid temperature of the gas is below some 
threshold. In the parameter file this is called :code:`temperature_threshold_K`
in which a good value is between 100 K and 1000 K. 

The fourth criterion is the maximal allowed density for gas. In the case that
we do not want gas to reach densities above a certain value we can 
instantaneously convert gas particles to stars if necessary, this variable is 
called :code:`threshold_max_density_H_p_cm3`, typical values for EAGLE like
models are around :math:`10^5` particles per cubic centimeter, in COLIBRE we 
therefore might need to use a higher value, but for now we are not using this
variable. This variable is optional and when not specified the value is set to
infinity.

The fifth and last criterion is the criterion that makes all gas above a certain
subgrid density have a star formation rate. Like the previous criterion, this 
criterion is optional and so when not specified is set to infinity. This 
criterion can be specified by using :code:`subgrid_density_threshold_H_p_CM3`.
In combination with the temperature model a value of :math:`10^2` is a good 
value. 

Lastly it is important to realize for models that use the subgrid temperature 
criteria that there is an implicit dependence with this criteria in the point
at which the cooling calculates the subgrid temperature. This is specified by
:code:`delta_logTEOS_subgrid_properties` which specifies how far away from the
EoS we calculate the subgrid temperature. A default value for this is 0.2, but 
using 0.3 or 0.5 makes no practical difference. 

+----------------------------------------+---------------------------------------+-----------------------+
| Name                                   | Description                           | Comments              |
+========================================+=======================================+=======================+
|| ``min_over_density``                  | | Is gas in a high enough over density| | Mandatory           |
+----------------------------------------+---------------------------------------+-----------------------+
|| ``alpha_virial``                      | | Is the virial criterion satisfied?  | | Mandatory           |
+----------------------------------------+---------------------------------------+-----------------------+
|| ``temperature_threshold``             | | Particles need to have a lower      | | Mandatory           |
|                                        | | subgrid temperature than this to be |                       |
|                                        | | star forming.                       |                       |
+----------------------------------------+---------------------------------------+-----------------------+
|| ``subgrid_density_threshold_H_p_CM3`` | | Subgrid density threshold from      | | Not specifying the  |
|                                        | | which to always be star forming     | | variable sets it    |
|                                        |                                       | | equal to infinity   |
+----------------------------------------+---------------------------------------+-----------------------+
|| ``threshold_max_density_H_p_cm3``     | | Maximal gas density, gas above this | | Not specifying the  |
|                                        | | density is instantaneously converted| | variable sets it    |
|                                        | | into stars.                         | | to infinity         |
+----------------------------------------+---------------------------------------+-----------------------+

Schmidt law
^^^^^^^^^^^

The first model that sets the star formation rate is the Schmidt law. This 
model can be specified by using :code:`SF_model:SchmidtLaw`. 

The star formation rate for the Schmidt law is implemented as:
:math:`\dot{\rho}_\star = \epsilon_{sfe} \frac{\rho}{\tau_{ff}}`. Which is 
implemented as :math:`\dot{m}_\star = \epsilon_{sfe} \frac{m_{gas}}{\tau_{ff}}`,
with the free fall time calculated as 
:math:`\tau_{ff} = \sqrt{\frac{3\pi}{32 G \rho}}`.
In which :math:`G` is the gravitational constant and :math:`\rho` is the SPH 
density. 

Given a star formation rate the probability of a converting the particle to a 
star is given by: :math:`Prob=\min\left(\frac{\dot{m}_*\Delta t}{m_g},1\right)`.

It seems clear that the Schmidt law only has one additional free parameter, the
star formation efficiency: :code:`star_formation_efficiency`, a good value for
the star formation efficiency is :math:`0.01`. 

Pressure law
^^^^^^^^^^^^

The second model that is available to use is the pressure law model from 
Schaye & Dalla Vecchia (2008) that was used in the EAGLE simulations. This
model sets the star formation rate based on input values for the 
Kennicutt-Schmidt relation. The star formation rate in this model is given
by:
:math:`\dot{\rho}_\star = A (1 M_\odot pc^{-2})^{-n} \left( \frac{\gamma}{G} f_g P_{tot} \right)^{(n-1)/2} \rho_g`

In this model we have 3 variables that need to be set. The first is the 
observed slope of the Kennicutt-Schmidt relation, this is specified in the 
parameter file as :code:`KS_exponent`, a typical value for this is 1.4.
The second is the normalization of the Kennicutt-Schmidt relation, which can
be set by :code:`KS_normalisation_Msun_p_yr_p_kpc2` and which typical has a 
value of :math:`1.515e-4`. The last variable is optional and is the gas 
fraction by default we assume the gas fraction is 1, but it can also be 
specified as an other value it can be specified by :code:`gas_fraction`.

Examples
^^^^^^^^

A run with all the parameters to run with a Schmidt law and a virial criterion
are given by:

.. code:: YAML

    # COLIBRE star formation parameters
    COLIBREStarFormation:
      min_over_density:                  57.7
      alpha_virial:                      1.
      temperature_threshold:             0
      SF_model:                          SchmidtLaw
      star_formation_efficiency:         0.01 

A run that is run with the pressure law and uses a non unity gas fraction, a 
temperature criterion and both density criteria looks like:

.. code:: YAML

    # COLIBRE star formation parameters
    COLIBREStarFormation:
      min_over_density:                  57.7
      alpha_virial:                      1.
      temperature_threshold:             1000
      threshold_max_density_H_p_cm3:     1e5
      subgrid_density_threshold_H_p_CM3: 1e2
      SF_model:                          PressureLaw
      KS_exponent:                       1.4
      KS_normalisation_Msun_p_yr_p_kpc2: 1.515e-4
      gas_fraction:                      0.3


.. _COLIBRE_delay_time_distributions:

SNIa Delay Time Distributions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the COLIBRE feedback there are several different DTDs implemented:

* The exponential DTD: used in the EAGLE simulations: 
  :math:`\frac{\nu}{\tau}~ \exp(-\frac{t_{init}}{\tau}) \exp(-\frac{t}{\tau})`
* A general power law: :math:`\nu ~\text{norm}~ t^{-\beta}`
* A power law with :math:`\beta=1`, :math:`\nu~ \text{norm}~ t^{-1}` 
* A constant plus Gaussian DTD, :math:`\frac{\nu_1}{t_{norm}} +
  \nu_2 ~ \frac{1}{\sqrt{2\pi \sigma^2}} \exp\left(-\frac{(t-t_{center})^2}{sigma^2}\right)`
* A constant DTD, :math:`\frac{\nu}{t_{norm}}`
* A broken power law DTD: :math:`\nu~ \text{norm}~ (t/t_b)^{-p}`, in which 
  :math:`p=p_{short}` for :math:`t<t_b` and :math:`p=p_{long}` for 
  :math:`t>t_b`.

A summary of the different DTDs is shown below which shows the 
free parameters in the model. 

+-----------------------+------------------------------------------------+
| Name                  | free parameters                                | 
+-----------------------+------------------------------------------------+
| | EAGLE or exponential| | SNIa_delay_time_Gyr, SNIa_efficiency_p_Msun, |
|                       | | SNIa_timescale_Gyr                           |
+-----------------------+------------------------------------------------+
| | power-law           | | SNIa_delay_time_Gyr, SNIa_efficiency_p_Msun, |
|                       | | normalization_timescale_Gyr, power_law_slope |
+-----------------------+------------------------------------------------+
| | power-law-beta-one  | | SNIa_delay_time_Gyr, SNIa_efficiency_p_Msun, |
|                       | | normalization_timescale_Gyr                  |
+-----------------------+------------------------------------------------+
| | gaussian            | | SNIa_delay_time_Gyr,                         |
|                       | | SNIa_efficiency_const_p_Msun,                |
|                       | | SNIa_efficiency_gauss_p_Msun,                |
|                       | | normalization_timescale_Gyr,                 |
|                       | | characteristic_time_Gyr,                     |
|                       | | STD_characteristic_time_Gyr                  |
+-----------------------+------------------------------------------------+
| | constant            | | SNIa_delay_time_Gyr, SNIa_efficiency_p_Msun, |
|                       | | normalization_timescale_Gyr                  |
+-----------------------+------------------------------------------------+
| | broken-power-law    | | SNIa_delay_time_Gyr, SNIa_efficiency_p_Msun, |
|                       | | power_law_slope_short_time,                  |
|                       | | power_law_slope_long_time,                   |
|                       | | break_time_Gyr, normalization_timescale_Gyr  |
+-----------------------+------------------------------------------------+

Example of how to run the code with the different DTDs:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In this part we show how to run with different DTDs in the COLIBRE 
subgrid model. The code by default runs with a power law with slope 
:math:`\beta=1`, the parameter file for the DTD in this case looks like:

.. code:: YAML
      
    #DTD parameters
    SNIaDTD:
      SNIa_efficiency_p_Msun:       0.001
      SNIa_delay_time_Gyr:          0.04
      normalization_timescale_Gyr:  13.6

Running the code with the exponential DTD is also possible as:

.. code:: YAML
      
    #DTD parameters
    SNIaDTD:
      SNIa_efficiency_p_Msun:       0.001
      SNIa_delay_time_Gyr:          0.04
      SNIa_timescale_Gyr:           2.0

Running the code with the general power law (:math:`\beta \neq 1`):

.. code:: YAML
      
    #DTD parameters
    SNIaDTD:
      SNIa_efficiency_p_Msun:       0.001
      SNIa_delay_time_Gyr:          0.04
      power_law_slope:              1.1
      normalization_timescale_Gyr:  13.6

Running the code with the constant + gaussian model:

.. code:: YAML
      
    #DTD parameters
    SNIaDTD:
      SNIa_delay_time_Gyr:          0.04
      normalization_timescale_Gyr:  13.6
      SNIa_efficiency_const_p_Msun: 0.001
      SNIa_efficiency_gauss_p_Msun: 0.001
      characteristic_time_Gyr:      4.0
      STD_characteristic_time_Gyr:  2.0
    
Running the code with the constant model:

.. code:: YAML
      
    #DTD parameters
    SNIaDTD:
      SNIa_efficiency_p_Msun:       0.001
      SNIa_delay_time_Gyr:          0.04
      normalization_timescale_Gyr:  13.6

Running the code with the broken power law:

.. code:: YAML
      
    #DTD parameters
    SNIaDTD:
      SNIa_efficiency_p_Msun:       0.001
      SNIa_delay_time_Gyr:          0.04
      normalization_timescale_Gyr:  13.6
      power_law_slope_short_time:   0.5
      power_law_slope_long_time:    1.1
      break_time_Gyr:               0.4

