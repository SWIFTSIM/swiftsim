.. CHIMES parameters 
   Alexander Richings 28th January 2020 

.. _CHIMES_param:

CHIMES parameters
-----------------

The following parameters are specified in the parameter file when running SWIFT with the CHIMES module: 

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
