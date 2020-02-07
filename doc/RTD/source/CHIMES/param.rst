.. CHIMES parameters 
   Alexander Richings 28th January 2020 

.. _CHIMES_param:

CHIMES parameters
-----------------

The parameters needed when running SWIFT with the CHIMES module are described in detail below. At the bottom of this page we have included example sets of CHIMES parameters to familiarise new users with typical parameter values for common set ups. 

+------------------------------------+---------------------------------------------------------------+
| Parameter                          | Description                                                   |
+====================================+===============================================================+
| ``data_path``                      | | Path to the chimes-data repository, containing the CHIMES   |
|                                    | | data files.                                                 |
|                                    |                                                               |
+------------------------------------+---------------------------------------------------------------+
| ``EqmAbundanceTable``              | | Path to the equilibrium abundance table, used when evolving |
|                                    | | the temperature with equilibrium cooling in CHIMES,         |
|                                    | | relative to the ``data_path`` path. If using  a redshift-   |
|                                    | | dependent UV background, this should be the path to the     |
|                                    | | directory containing the redshift-dependent equilibrium     |
|                                    | | tables. Otherwise, it should point to the HDF5 table        |
|                                    | | itself.                                                     |
|                                    |                                                               |
+------------------------------------+---------------------------------------------------------------+
| ``UV_field_flag``                  | | Integer flag that sets the UV radiation field to be used in |
|                                    | | the CHIMES chemistry solver. Possible values are as         |
|                                    | | follows:                                                    |
|                                    | | 0 - No UV radiation.                                        |
|                                    | | 1 - Single user-defined UV spectrum (see below).            |
|                                    | | 2 - Two UV spectra, corresponding to the UV background plus |
|                                    | |     a local ISRF, scaled according to Ploeckinger et al.    |
|                                    | |     (in prep).                                              |
|                                    |                                                               |
+------------------------------------+---------------------------------------------------------------+
| ``use_redshift_dependent_UVB``     | | Integer flag to determine whether to interpolate the UV     |
|                                    | | background to the current redshift. NOTE if the             |
|                                    | | ``UV_field_flag`` parameter has been set to zero, this      |
|                                    | | parameter must also be set to zero, as no UVB is included.  |
|                                    | | Possible values are as follows:                             |
|                                    | | 0 - No interpolation with redshift. All UV radiation fields |
|                                    | |     (if included) are fixed.                                |
|                                    | | 1 - Interpolate the photoionisation cross sections tables   |
|                                    | |     and the strength of the UVB to the current redshift. If |
|                                    | |     ``UV_field_flag == 2``, it is only the UVB component    |
|                                    | |     that is interpolated in this way.                       |
|                                    | | 2 - Interpolate both the UVB cross sections tables and the  |
|                                    | |     equilibrium abundance tables to the current redshift.   |
|                                    |                                                               |
+------------------------------------+---------------------------------------------------------------+
| ``rad_field_norm_factor``          | | Multiplicative factor to scale the normalisation of the UV  |
|                                    | | radiation field. The behaviour of this parameter depends on |
|                                    | | the ``UV_field_flag`` as follows:                           |
|                                    | | ``UV_field_flag == 0`` - This parameter is not used.        |
|                                    | | ``UV_field_flag == 1`` - Scales the single radiation field  |
|                                    | |    specified by the user, relative to the default strength  |
|                                    | |    of that field in the cross sections table (see below).   |
|                                    | | ``UV_field_flag == 2`` - Scales only the normalisation of   |
|                                    | |    the local ISRF component, relative to a default Milky    |
|                                    | |    Way scaling. Note that the default COLIBRE model uses a  |
|                                    | |    value of ``0.1``, corresponding to reducing the Milky    |  
|                                    | |    Way ISRF by a factor of 10 (see Ploeckinger et al.       | 
|                                    | |    in prep).                                                |
|                                    |                                                               |
+------------------------------------+---------------------------------------------------------------+
| ``PhotoIonTable``                  | | Only used when ``UV_field_flag == 1``.                      |
|                                    | | Specifies the photoionisation cross sections table for the  |
|                                    | | single UV radiation field specified by the user. If using a |
|                                    | | redshift-dependent UVB, this parameter gives the path to    |
|                                    | | the directory containing the UVB cross sections tables.     |
|                                    | | Otherwise, this parameter gives the path to the HDF5 table  |
|                                    | | itself. All paths are relative to the ``chimes_data`` path. |
|                                    |                                                               |
+------------------------------------+---------------------------------------------------------------+
| ``PhotoIonTable_UVB``              | | Only used when ``UV_field_flag == 2``.                      |
|                                    | | Specifies the photoionisation cross sections table for the  |
|                                    | | UVB component of the radiation field. If using a redshift-  |
|                                    | | dependent UVB, this parameter gives the path to the         |
|                                    | | directory containing the UVB cross sections tables.         |
|                                    | | Otherwise, this parameter gives the path to the HDF5 table  |
|                                    | | itself. All paths are relative to the ``chimes_data`` path. |
|                                    |                                                               |
+------------------------------------+---------------------------------------------------------------+
| ``PhotoIonTable_ISRF``             | | Only used when ``UV_field_flag == 2``.                      |
|                                    | | Specifies the photoionisation cross sections table for the  |
|                                    | | local ISRF component of the radiation field. This parameter |
|                                    | | gives the path to the HDF5 table itself, relative to the    |
|                                    | | ``chimes_data`` path.                                       |
|                                    |                                                               |
+------------------------------------+---------------------------------------------------------------+
| ``Shielding_flag``                 | | Integer flag that sets the local shielding approximation.   |
|                                    | | Possible values are as follows:                             |
|                                    | | 0 - No shielding.                                           |
|                                    | | 1 - Include shielding with a shielding length equal to the  |
|                                    | |     local Jeans length.                                     |
|                                    | | 2 - Include shielding with a shielding length scaled        |
|                                    | |     according to Ploeckinger et al. (in prep), as in        |
|                                    | |     COLIBRE.                                                |
|                                    |                                                               |
+------------------------------------+---------------------------------------------------------------+
| ``shielding_length_factor``        | | Only used if ``Shielding_flag > 0`` or                      |
|                                    | | ``UV_field_flag == 2``.                                     |
|                                    | | Multiplicative factor to scale the local shielding length.  |
|                                    | | If ``UV_field_flag == 2``, this is also used to scale the   |
|                                    | | reference column density, N_ref, for calculating the local  |
|                                    | | ISRF (see Ploeckinger et al. in prep).                      |
|                                    |                                                               |
+------------------------------------+---------------------------------------------------------------+
| ``max_shielding_length``           | | Only used if ``Shielding_flag > 0`` or                      |
|                                    | | ``UV_field_flag == 2``.                                     |
|                                    | | Maximum shielding length, in code units. If                 |
|                                    | | ``UV_field_flag == 2``, this is also used to limit the      |
|                                    | | shielding length used for the reference column density,     |
|                                    | | N_ref, in calculating the local ISRF (see Ploeckinger et    |
|                                    | | al. in prep).                                               |
|                                    |                                                               |
+------------------------------------+---------------------------------------------------------------+
| ``H_reion_z``                      | | Redshift of reionisation.                                   |
|                                    |                                                               |
+------------------------------------+---------------------------------------------------------------+
| ``init_abundance_mode``            | | Integer flag that determines how the initial CHIMES         |
|                                    | | abundances are set at the beginning of the simulation.      |
|                                    | | Possible values are as follows:                             |
|                                    | | 0 - Set each element to one ionisation state (e.g. neutral, |
|                                    | |     single ionised etc.) according to the ``InitIonState``  |
|                                    | |     parameter (see below).                                  |
|                                    | | 1 - Read in the abundances for each particle from           |
|                                    | |     pre-computed equilibrium abundance tables, as a         |
|                                    | |     function of temperature, density and metallicity.       |
|                                    | | 2 - Compute the initial equilibrium abundances for each     |
|                                    | |     particle at the beginning of the simulation.            |
|                                    |                                                               |
+------------------------------------+---------------------------------------------------------------+
| ``InitIonState``                   | | Only used if ``init_abundance_mode == 0``.                  |
|                                    | | Gives the initial ionisation state to set each element to   |
|                                    | | at the beginning of the simulation (0: neutral, 1: single   |
|                                    | | ionised etc.).                                              |
|                                    |                                                               |
+------------------------------------+---------------------------------------------------------------+
| ``cosmic_ray_rate``                | | Ionisation rate of HI by cosmic rays, in units of s^-1.     |
|                                    | | If ``UV_field_flag == 2``, this needs to be set to the      |
|                                    | | Milky Way value, i.e. ``1.8e-16``. It will then be scaled   |
|                                    | | along with the local ISRF as in Ploeckinger et al. (in      |
|                                    | | prep). Note that, in this case, the Milky Way value is      |
|                                    | | then multiplied by ``rad_field_norm_factor``. So if you are |
|                                    | | using the fiducial cooling model from Ploeckinger et al (in |
|                                    | | prep), i.e. with ``rad_field_norm_factor = 0.1``, you do    |
|                                    | | not need to reduce ``cosmic_ray_rate`` by another factor of |
|                                    | | 10, as this is already dealt with internally in the code.   |
|                                    |                                                               |
+------------------------------------+---------------------------------------------------------------+
| ``relativeTolerance``              | | Relative tolerance parameter, used to define the accuracy   |
|                                    | | of the chemistry and cooling integration in CHIMES.         |
|                                    |                                                               |
+------------------------------------+---------------------------------------------------------------+
| ``absoluteTolerance``              | | Absolute tolerance parameter, used to define the accuracy   |
|                                    | | of the chemistry integration in CHIMES.                     |
|                                    |                                                               |
+------------------------------------+---------------------------------------------------------------+
| ``thermalAbsoluteTolerance``       | | Absolute tolerance parameter for the temperature, used to   |
|                                    | | define the accuracy of the cooling integration in CHIMES.   |
|                                    |                                                               |
+------------------------------------+---------------------------------------------------------------+
| ``explicitTolerance``              | | Tolerance parameter that determines when we just use the    |
|                                    | | explicit solution for the chemistry and cooling integration |
|                                    | | in CHIMES.                                                  |
|                                    |                                                               |
+------------------------------------+---------------------------------------------------------------+
| ``scale_metal_tolerances``         | | Integer flag (``0`` or ``1``) that determines whether to    |
|                                    | | scale the abolute tolerances of each species by its         |
|                                    | | corresponding element abundance.                            |
|                                    |                                                               |
+------------------------------------+---------------------------------------------------------------+
| ``T_mol``                          | | Maximum temperature for the molecular network. Above this   |
|                                    | | temperature we skip over all reactions involving molecules  |
|                                    | | and set all molecule abundances to zero.                    |
|                                    |                                                               |
+------------------------------------+---------------------------------------------------------------+
| ``ChemistryEqmMode``               | | Integer flag to set the equilibrium mode in CHIMES.         |
|                                    | | Possible values are as follows:                             |
|                                    | | 0 - Evolve CHIMES species in non-equilibrium.               |
|                                    | | 1 - Evolve the cooling with CHIMES species in equilibrium,  |
|                                    | |     i.e. use the pre-computed equilibrium abundance tables. |
|                                    |                                                               |
+------------------------------------+---------------------------------------------------------------+
| ``ThermEvolOn``                    | | Integer flag to switch the temperature evolution on/off in  |
|                                    | | CHIMES. Possible values are as follows:                     |
|                                    | | 0 - Evolve the chemical abundances at fixed temperature.    |
|                                    | | 1 - Evolve both the temperature and the chemical            |
|                                    | |     abundances.                                             |
|                                    |                                                               |
+------------------------------------+---------------------------------------------------------------+
| ``chimes_debug``                   | | Integer flag for additional debug output from CHIMES.       |
|                                    | | Possible values are as follows:                             |
|                                    | | 0 - No additional output.                                   |
|                                    | | 1 - If CVode returns a non-zero flag (i.e. it returns a     |
|                                    | |     CVode error or warning), print out all of the variables |
|                                    | |     in the ChimesGasVars stucture.                          |
|                                    |                                                               |
+------------------------------------+---------------------------------------------------------------+
| ``S_over_Si_in_solar``             | | S to Si ratio relative to the Solar ratio.                  |
|                                    |                                                               |
+------------------------------------+---------------------------------------------------------------+
| ``Ca_over_Si_in_solar``            | | Ca to Si ratio relative to the Solar ratio.                 |
|                                    |                                                               |
+------------------------------------+---------------------------------------------------------------+
| ``colibre_metal_depletion``        | | Integer flag (``0`` or ``1``) that sets whether to reduce   |
|                                    | | the gas-phase element abundances in CHIMES due to dust      | 
|                                    | | depletion according to Ploeckinger et al. (in prep).        |
|                                    |                                                               |
+------------------------------------+---------------------------------------------------------------+
|``delta_logTEOS_subgrid_properties``| | Distance from the EOS to use the thermal equilibrium        |
|                                    | | temperature for subgrid properties, and to evolve the       |
|                                    | | cooling using equilibrium abundances.                       |
|                                    |                                                               |
+------------------------------------+---------------------------------------------------------------+
| ``use_colibre_subgrid_EOS``        | | Integer flag (``0`` or ``1``) to use the subgrid density    |
|                                    | | and temperature from the COLIBRE cooling tables (see        |
|                                    | | Ploeckinger et al. in prep).                                |
|                                    |                                                               |
+------------------------------------+---------------------------------------------------------------+
| ``rapid_cooling_threshold``        | | Threshold in ``dt / t_cool`` above which we are in the      |
|                                    | | rapid cooling regime (i.e. cooling time is short compared   |
|                                    | | to the hydro time-step). In this case, the cooling routines |
|                                    | | will instantly update the particle's internal energy with   |
|                                    | | the final internal energy at the end of the time-step.      |
|                                    | | If the particle is below this threshold, its internal       |
|                                    | | energy will instead be drifted using the ``du_dt`` from the |
|                                    | | cooling routines.                                           |
|                                    | | If this parameter is negative, we always drift the          |
|                                    | | temperature.                                                |
|                                    |                                                               |
+------------------------------------+---------------------------------------------------------------+
| ``HIIregion_ionization_fraction``  | | Ionisation fraction of HII regions (``0`` - fully neutral;  |
|                                    | | ``1`` - fully ionised).                                     |
|                                    |                                                               |
+------------------------------------+---------------------------------------------------------------+
| ``HIIregion_temperature``          | | Minimum temperature of HII regions.                         |
|                                    |                                                               |
+------------------------------------+---------------------------------------------------------------+
| ``use_hybrid_cooling``             | | Integer flag to set whether to use the hybrid cooling mode. |
|                                    | | Possible values are as follows:                             |
|                                    | | 0 - Don't use hybrid cooling. If any elements are switched  |
|                                    | |     off in the CHIMES network, their cooling is neglected.  |
|                                    | | 1 - Use hybrid cooling. If any elements are switched off in |
|                                    | |     the CHIMES network, we look up their cooling rates from |
|                                    | |     the COLIBRE cooling tables and add it on to the total   |
|                                    | |     cooling rate from CHIMES when integrating the           |
|                                    | |     temperature.                                            |
|                                    |                                                               |
+------------------------------------+---------------------------------------------------------------+
| ``IncludeCarbon``                  | | Integer flag (``0`` or ``1``) to set whether to include     |
|                                    | | carbon in the CHIMES network.                               |
|                                    |                                                               |
+------------------------------------+---------------------------------------------------------------+
| ``IncludeNitrogen``                | | Integer flag (``0`` or ``1``) to set whether to include     |
|                                    | | nitrogen in the CHIMES network.                             |
|                                    |                                                               |
+------------------------------------+---------------------------------------------------------------+
| ``IncludeOxygen``                  | | Integer flag (``0`` or ``1``) to set whether to include     |
|                                    | | oxygen in the CHIMES network.                               |
|                                    |                                                               |
+------------------------------------+---------------------------------------------------------------+
| ``IncludeNeon``                    | | Integer flag (``0`` or ``1``) to set whether to include     |
|                                    | | neon in the CHIMES network.                                 |
|                                    |                                                               |
+------------------------------------+---------------------------------------------------------------+
| ``IncludeMagnesium``               | | Integer flag (``0`` or ``1``) to set whether to include     |
|                                    | | magnesium in the CHIMES network.                            |
|                                    |                                                               |
+------------------------------------+---------------------------------------------------------------+
| ``IncludeSilicon``                 | | Integer flag (``0`` or ``1``) to set whether to include     |
|                                    | | silicon in the CHIMES network.                              |
|                                    |                                                               |
+------------------------------------+---------------------------------------------------------------+
| ``IncludeSulphur``                 | | Integer flag (``0`` or ``1``) to set whether to include     |
|                                    | | sulphur in the CHIMES network.                              |
|                                    |                                                               |
+------------------------------------+---------------------------------------------------------------+
| ``IncludeCalcium``                 | | Integer flag (``0`` or ``1``) to set whether to include     |
|                                    | | calcium in the CHIMES network.                              |
|                                    |                                                               |
+------------------------------------+---------------------------------------------------------------+
| ``IncludeIron``                    | | Integer flag (``0`` or ``1``) to set whether to include     |
|                                    | | iron in the CHIMES network.                                 |
|                                    |                                                               |
+------------------------------------+---------------------------------------------------------------+
| ``colibre_table_path``             | | Only used if ``use_hybrid_cooling == 1`` or                 |
|                                    | | ``use_colibre_subgrid_EOS == 1``. Specifies the path to the |
|                                    | | COLIBRE cooling table (see Ploeckinger et al. in prep).     |
|                                    |                                                               |
+------------------------------------+---------------------------------------------------------------+


Example: Isolated galaxy hybrid cooling
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following set of CHIMES parameters are suitable for running a non-cosmological isolated galaxy with hybrid cooling, i.e. with hydrogen and helium in non-equilibrium and metal cooling read in from the COLIBRE cooling tables in equilibrium. This example uses the COLIBRE models for the UV radiation field, shielding and metal depletion. Note that you will need to change the `data_path` and `colibre_table_path` parameters to point to the correct directory on your system. 

.. code:: YAML

    # CHIMES cooling parameters
    CHIMESCooling: 
      data_path:                  /path/to/chimes-data 
      EqmAbundanceTable:          colibre_HHe/z0.000_eqm.hdf5 
      PhotoIonTable_UVB:          HM12_cross_sections/z0.000_cross_sections.hdf5 
      PhotoIonTable_ISRF:         cross_sections_B87.hdf5 
      UV_field_flag:              2 
      Shielding_flag:             2
      use_redshift_dependent_UVB: 0
      shielding_length_factor:    0.5 
      max_shielding_length:       100.0 
      rad_field_norm_factor:      0.1 
      init_abundance_mode:        1 
      colibre_metal_depletion:    1
      relativeTolerance:          1e-4 
      absoluteTolerance:          1e-10 
      thermalAbsoluteTolerance:   1e-40 
      explicitTolerance:          0.05 
      scale_metal_tolerances:     1 
      T_mol:                      1.0e5 
      ChemistryEqmMode:           0 
      ThermEvolOn:                1 
      chimes_debug:               0 
      cosmic_ray_rate:            1.8e-16 
      delta_logTEOS_subgrid_properties:  0.2 
      use_colibre_subgrid_EOS:    1 
      use_hybrid_cooling:         1 
      rapid_cooling_threshold:    1.0 
      HIIregion_ionization_fraction: 1.0 
      HIIregion_temperature:      1.0e4 
      colibre_table_path:         /path/to/UV_dust1_CR1_G1_shield1.hdf5 
      H_reion_z:                  7.5 
      S_over_Si_in_solar:         1.0
      Ca_over_Si_in_solar:        1.0
      IncludeCarbon:              0 
      IncludeNitrogen:            0 
      IncludeOxygen:              0 
      IncludeNeon:                0 
      IncludeMagnesium:           0 
      IncludeSilicon:             0 
      IncludeSulphur:             0 
      IncludeCalcium:             0 
      IncludeIron:                0 


Example: Cosmological box hybrid cooling
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The following set of CHIMES parameters are suitable for running a cosmological box with hybrid cooling, i.e. with hydrogen and helium in non-equilibrium and metal cooling read in from the COLIBRE cooling tables in equilibrium. This example uses the COLIBRE models for the UV radiation field, shielding and metal depletion. Note that you will need to change the `data_path` and `colibre_table_path` parameters to point to the correct directory on your system. 

Compared to the isolated galaxy example above, this example uses a redshift-dependent UVB. 

.. code:: YAML

    # CHIMES cooling parameters
    CHIMESCooling: 
      data_path:                  /path/to/chimes-data 
      EqmAbundanceTable:          colibre_HHe 
      PhotoIonTable_UVB:          HM12_cross_sections
      PhotoIonTable_ISRF:         cross_sections_B87.hdf5 
      UV_field_flag:              2 
      Shielding_flag:             2
      use_redshift_dependent_UVB: 2
      shielding_length_factor:    0.5 
      max_shielding_length:       100.0 
      rad_field_norm_factor:      0.1 
      init_abundance_mode:        1 
      colibre_metal_depletion:    1
      relativeTolerance:          1e-4 
      absoluteTolerance:          1e-10 
      thermalAbsoluteTolerance:   1e-40 
      explicitTolerance:          0.05 
      scale_metal_tolerances:     1 
      T_mol:                      1.0e5 
      ChemistryEqmMode:           0 
      ThermEvolOn:                1 
      chimes_debug:               0 
      cosmic_ray_rate:            1.8e-16 
      delta_logTEOS_subgrid_properties:  0.2 
      use_colibre_subgrid_EOS:    1 
      use_hybrid_cooling:         1 
      rapid_cooling_threshold:    1.0 
      HIIregion_ionization_fraction: 1.0 
      HIIregion_temperature:      1.0e4 
      colibre_table_path:         /path/to/UV_dust1_CR1_G1_shield1.hdf5 
      H_reion_z:                  7.5 
      S_over_Si_in_solar:         1.0
      Ca_over_Si_in_solar:        1.0
      IncludeCarbon:              0 
      IncludeNitrogen:            0 
      IncludeOxygen:              0 
      IncludeNeon:                0 
      IncludeMagnesium:           0 
      IncludeSilicon:             0 
      IncludeSulphur:             0 
      IncludeCalcium:             0 
      IncludeIron:                0 
