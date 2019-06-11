.. COLIBRE sub-grid model
   Matthieu Schaller, 20th December 2018
   Folkert Nobels, 3th of June


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

The COLIBRE star formation has optional criteria for the gas to be 
star forming a list of optional paramters that can be set are given by:

+---------------------------------------+---------------------------------------+-----------------------+
| Name                                  | Description                           | Comments              |
+=======================================+=======================================+=======================+
| | ``temperature_threshold``           | | Particles need to have a lower      | | Not specifying the  |
|                                       | | temerpature than this to be         | | variables sets it   |
|                                       | | starforming.                        | | equal to infinity   |
+---------------------------------------+---------------------------------------+-----------------------+
| | ``EOS_temperature_margin_dex``       | | Off set in dex from the equation of | | Not specifying the  |
|                                       | | state for gas to form stars         | | variable sets it    |
|                                       |                                       | | equal to 10.        | 
+---------------------------------------+---------------------------------------+-----------------------+ 
| | ``threshold_max_density_H_p_cm3``   | | Maximal gas density, gas above this | | Not specifying the  |
|                                       | | density is instantaneously converted| | variable sets it    |
|                                       | | into stars.                         | | to infinity         |
+---------------------------------------+---------------------------------------+-----------------------+ 

A run with all the paramters will have a YAML file that looks like:

.. code:: YAML

    # COLIBRE star formation parameters
    COLIBREStarFormation:
      min_over_density:                 57.7
      star_formation_efficiency:        0.01 
      temperature_threshold:            1000
      threshold_max_density_H_p_cm3:    1e5
      EOS_temperature_margin_dex:        0.2

Code that only has a temperature margin dex critiria will look like:

.. code:: YAML

    # COLIBRE star formation parameters
    COLIBREStarFormation:
      min_over_density:                 57.7
      star_formation_efficiency:        0.01 
      threshold_max_density_H_p_cm3:    1e5
      EOS_temperature_margin_dex:        0.2

And code with only a temperature threshold will look like:

.. code:: YAML

    # COLIBRE star formation parameters
    COLIBREStarFormation:
      min_over_density:                 57.7
      star_formation_efficiency:        0.01 
      temperature_threshold:            1000
      threshold_max_density_H_p_cm3:    1e5

In the future new additional star formation criteria will be added to this 
routine.
