.. CHIMES config 
   Alexander Richings 28th January 2020 

.. _CHIMES_config:

Configuring SWIFT with CHIMES
------------------------------

To build SWIFT with the CHIMES module, you will need the following configuration options: 

+------------------------------------+---------------------------------------------------------------+
| Configuration option               | Description                                                   |
+====================================+===============================================================+
| ``--with-cooling=CHIMES``          | | Sets the SWIFT cooling function to use the CHIMES module.   |
|                                    |                                                               |
+------------------------------------+---------------------------------------------------------------+
| ``--with-sundials=/sundials/path`` | | Path to the directory where the Sundials library is         |
|                                    | | installed.                                                  |
|                                    |                                                               |
+------------------------------------+---------------------------------------------------------------+
| ``--with-chimes-network-size=157`` | | The number of ions and molecules included in the CHIMES     |
|                                    | | network, i.e. that are evolved in non-equilibrium. This has |
|                                    | | to match the size of the network determined from the        |
|                                    | | elements that are included according to the                 |
|                                    | | Include[Element] parameters (see the parameters section).   |
|                                    | | Common values:                                              |
|                                    | | H + He only (as used in COLIBRE): 10                        |
|                                    | | Full network: 157                                           |
|                                    |                                                               |
+------------------------------------+---------------------------------------------------------------+
