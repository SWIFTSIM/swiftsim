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
| ``--with-sundials=/sundials/path`` | | Builds SWIFT with the Sundials libraries, which are         |
|                                    | | required by the CHIMES module. If you just use              |
|                                    | | ``--with-sundials``, it will try to find the Sundials       |
|                                    | | libraries automatically on your system. However, if it      |
|                                    | | cannot find these libraries automatically (for example, if  |
|                                    | | you built the libraries yourself), then you can directly    |
|                                    | | pass it the path to the directory where the Sundials        |
|                                    | | libraries are installed.                                    |
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


