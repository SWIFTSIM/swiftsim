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

The star formation is implemented as a simple density Schmidt law given
by :math:`\dot{\rho}_\star = \epsilon_{sfe} \frac{\rho}{\tau_{ff}}`. Which is 
implemented as :math:`\dot{m}_\star = \epsilon_{sfe} \frac{m_{gas}}{\tau_{ff}}`,
with the free fall time calculated as 
:math:`\tau_{ff} = \sqrt{\frac{3\pi}{32 G \rho}}`.
In which :math:`G` is the gravitational constant and :math:`\rho` is the SPH 
density. 

Given a star formation rate the probability of a converting the particle to a 
star is given by: :math:`Prob=\min\left(\frac{\dot{m}_*\Delta t}{m_g},1\right)`.

To prevent star formation in non-collapsed objects (for instance at high
redshift when the whole Universe has a density above the threshold), we apply an
over-density criterion. Only gas with a density larger than a multiple of the
critical density for closure can form stars.

The COLIBRE star formation has a few parameters that determine of the gas is 
star forming or not the list is shown below:

+----------------------------------------+---------------------------------------+-----------------------+
| Name                                   | Description                           | Comments              |
+========================================+=======================================+=======================+
|| ``temperature_threshold``             | | Particles need to have a lower      | | Mandatory           |
|                                        | | subgrid temperature than this to be |                       |
|                                        | | starforming.                        |                       |
+----------------------------------------+---------------------------------------+-----------------------+
|| ``subgrid_density_threshold_H_p_CM3`` | | Subgrid density threshold from      | | Not specifying the  |
|                                        | | which to always be starforming      | | variable sets it    |
|                                        |                                       | | equal to infinity   |
+----------------------------------------+---------------------------------------+-----------------------+
|| ``threshold_max_density_H_p_cm3``     | | Maximal gas density, gas above this | | Not specifying the  |
|                                        | | density is instantaneously converted| | variable sets it    |
|                                        | | into stars.                         | | to infinity         |
+----------------------------------------+---------------------------------------+-----------------------+

A run with all the paramters will have a YAML file that looks like:

.. code:: YAML

    # COLIBRE star formation parameters
    COLIBREStarFormation:
      min_over_density:                  57.7
      star_formation_efficiency:         0.01 
      temperature_threshold:             1000
      threshold_max_density_H_p_cm3:     1e5
      subgrid_density_threshold_H_p_CM3: 1e2

Code that only has a temperature threshold will look like:

.. code:: YAML

    # COLIBRE star formation parameters
    COLIBREStarFormation:
      min_over_density:                  57.7
      star_formation_efficiency:         0.01 
      temperature_threshold:             1000
      threshold_max_density_H_p_cm3:     1e5

In the future new additional star formation criteria will be added to this 
routine.


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

