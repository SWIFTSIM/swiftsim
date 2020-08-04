.. Lossy compression filters

.. _Compression_filters:

Compression Filters
~~~~~~~~~~~~~~~~~~~

Filters to compress the data in snapshots can be applied to reduce the
disk footprint of the datasets. The filters provided by SWIFT are
filters natively provided by HDF5, implying that the library will
automatically and transparently apply the reverse filter when reading
the data stored on disk. They can be applied in combination with, or
instead of, the lossless gzip compression filter.

**These compression filters are lossy.**

*This means that they will reduce the accuracy of the data stored. No check
is made inside SWIFT to verify that the applied filters make sense. Poor
choices can lead to all the values of a given array reduced to 0,
Inf, or to have lost too much accuracy to be useful. The onus is
entirely on the user to choose wisely how they want to compress the data.*

The available filters are listed below.

Scaling filters
---------------

The D-scale filters can be used to round floating-point values to a fixed
absolute accuracy.

They start by computing the minimum of an array that is then deducted from
all the values. The array is then multiplied by :math:`10^n` and truncated
to the nearest integer. These integers are stored with the minimal number
of bits required to store the values. That process is reversed when reading
the data.

For an array of values

+--------+--------+-------+
|  1.2345| -0.1267| 0.0897|
+--------+--------+-------+

and :math:`n=2`, we get stored on disk (but hidden to the user):

+--------+--------+-------+
|    136 |      0 |     22|
+--------+--------+-------+

This can be stored with 8 bits instead of the 32 bits needed to store the
original values in floating-point precision, realising a gain of 4x.

When reading the values (for example via ``h5py`` or ``swiftsimio``), that
process is transparently reversed and we get:

+--------+--------+-------+
|  1.2333| -0.1267| 0.0933|
+--------+--------+-------+

Using a scaling of :math:`n=2` hence rounds the values to two digits after
the decimal point.

SWIFT implements 4 variants of this filter:

 * ``DScale1`` --> Scale by :math:`10^1`
 * ``DScale2`` --> Scale by :math:`10^2`
 * ``DScale3`` --> Scale by :math:`10^3`
 * ``DScale6`` --> Scale by :math:`10^6`

An example application is to store the positions with `pc` accuracy in
simulations that use `Mpc` as their base unit by using the ``DScale6``
filter. 

