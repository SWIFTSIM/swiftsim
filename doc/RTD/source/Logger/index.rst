Logger Output
=============

The logger is a particle based output (e.g. snapshot) that takes into account the large difference of timescale.
If you have any question, a slack channel is available for it in SWIFT's slack.

To run it, you will need to use the configuration option ``--enable-logger``.
Currently the logger is implemented only for Gadget2 and the default gravity / stars, but can be easily extended to the other schemes by adding the logger structure to the particles (see ``src/hydro/Gadget2/hydro_part.h``).
The main parameters of the logger are ``Logger:delta_step`` and ``Logger:index_mem_frac`` that define the time accuracy of the logger and the number of index files.
The first parameter defines the number of active steps that a particle is doing before writing and the second defines the total storage size of the index files as function of the dump file.

Unfortunately, the API is not really developed yet. Therefore if you wish to dump another field, you will need to trick the logger by replacing a field in the ``logger_log_part`` function.

For reading, the python wrapper is available through the configuration option ``--with-python``. Once compiled, you will be able to use the file ``logger/examples/reader_example.py``.
The first argument is the basename of the index file and the second one is the time requested.
During the first reading, the library is manipulating the dump file and therefore it should not be killed and may take a bit more time than usual.


How to add a new field
----------------------

Only two files need modifications when adding a new field.
They define the behavior for the writer (e.g. ``src/gravity/MultiSoftening/gravity_logger.h``) and the reader
(e.g. ``logger/gravity/MultiSoftening/logger_gravity.h``).
In the two files, it is strongly recommended to always use the same order between the fields.
A few functions really need this order and some other do not, but I will not do the list.
Let's go over a simple example with the field ``new_field`` that is supposed to be already implemented in the ``gparts`` as a float.
We will call it ``NewField`` in the logger and in python.

First let's take a look at the writer.
When starting a simulation, the list of possible masks (value and data size) needs to be initialized:

.. code-block:: c

    INLINE static int gravity_logger_init(...) {
      /* Previous code */
      mask_data[n] = logger_add_field_to_logger("NewField", sizeof(float));
      return n + 1;
    }

Here n is the number of fields written before the addition of the new field.
The first step when writing a particle inside the logfile is to compute the size of the record.
At the same time, we can compute the record's mask:
 
.. code-block:: c

   INLINE static void gravity_logger_prepare_to_write_particle(...) {
     /* Previous code */
     *mask |= logger_add_field_to_mask(masks[1], "NewField", buffer_size);
   }

Finally, we need to write the particle according to the mask previously computed:

.. code-block:: c

   INLINE static char *gravity_logger_write_particle(...) {
     /* Previous code */

     /* First check that the field is supposed to be written */
     if (logger_should_write_field(mask_data[n], mask, "NewField")) {
       /* Copy the field to the buffer and update it */
       memcpy(buff, gp->new_field, sizeof(float));
       buff += sizeof(float);
     }
   }


Now we need to do the same for the reader.
First the field needs to be added to the reader's particle:

.. code-block:: c

   struct logger_gparticle {
     /* Previous fields */
     float new_field;
   }


Before reading the first particle, the code checks that the fields are the same (both in name and size)
in the logfile and the code. This is done in the following function:

.. code-block:: c

   __attribute__((always_inline)) INLINE static void logger_gparticle_check_fields(...) {
     /* Previous code */
     /* Check that the new field is present in the logfile's header */
     else if (strcmp(head->masks[i].name, "NewField") == 0) {
       /* Get the size of the field */
       size = sizeof(part.new_field);
     }
   }


If you wish, you can update the print function and initialize the new field in order
to have a safe behavior when the field is not read.
Then you will need to read the new field from the buffer:

.. code-block:: c

   __attribute__((always_inline)) INLINE static void logger_gparticle_read_field(...)
     /* Previous code */
     /* First check that the field is present in the buffer */
     } else if (strcmp("NewField", field) == 0) {
       /* Then copy the memory from the buffer */
       memcpy(&part->new_field, buff, size);
     }
   }


Once the particle is read, the code interpolates it to the required time.
Here you can use for example a linear interpolation or an Hermite interpolation if you have
at least one derivative.

.. code-block:: c

   __attribute__((always_inline)) INLINE static void logger_gparticle_interpolate(...) {

     /* Previous code */
     const float tmp = (part_next->new_field - part_curr->new_field);
     part_curr->new_field += tmp * scaling;
   }


Finally, if you wish to use python, you will need to expose it to python:

.. code-block:: c

   INLINE static int logger_gparticles_generate_python(...) {
     /* Previous code */
     /* f8 is the numpy string corresponding to a float. */
     list[n] = logger_loader_python_field("NewField", part, new_field, "f8");
     return n + 1;
   }
