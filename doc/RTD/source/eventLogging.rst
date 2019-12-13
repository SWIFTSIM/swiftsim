.. Event logger 
   Folkert Nobels, 9th of December 2019

Event logging in SWIFT
======================

In big simulations we only store a couple of snapshots during the whole 
simulation, some things require or need a higher time sampling. For this
there is event logging in SWIFT. Event logging takes a certain process
like SNII and on every node stores the amount of SNII events that took
place and combines them at preset times to write them to a log file. 
In SWIFT this can be done for any process as required/desired for the 
simulations. 

Adding a new event to log
~~~~~~~~~~~~~~~~~~~~~~~~~

To create a new event to log you needs to create a new file in the 
directory of the event logger, for example ``src/even_logger/event_logger_your_event.h``
in this file the first thing that you need to add is a ``event_history_your_event``
struct that contains a core ``event_history_logger`` struct that takes care
of all the logging essentials like when to log and writing to the file.
Besides this the ``event_history_your_struct`` should contain all the things you
actually want to log, like for example the number of events or the 
amount of mass or anything else. 

After constructing the struct we need to create the required functions,
the first function is the ``event_logger_your_event_init_log_file`` that is the 
function that initializes the log file with the used units and the different
columns in the log file. The second function is the ``event_logger_your_event_init``
function that initializes your event logger structure on each node and reads
the delta time to log. The third function is the ``event_logger_your_event_time_step``
which updates the event logging core after a time step. The fourth function
is the ``event_logger_your_event_log_data_general``, which is a function that 
actually logs the data to the log file. The ``functions event_logger_your_event_log_data``
and ``event_logger_your_event_log_data_end`` are both functions that log to the 
logger file, the first is called in general and the last only in the case of 
the end of the simulation and also closes the logging files. Further
``event_logger_your_event_log_event`` logs an event happening and locks the global
logging struct on your node such that not multiple threads access the variable,
and last ``event_logger_your_event_MPI_Reduce`` checks if we are on a logging step
and if so does a MPI_Reduce with all nodes to collect all the events on node 0
to log them to the log file.

If all these functions are written than the file should be included in 
event_logger.h and all the functions should be called in this file. 
Most of these have a very simular name as the ones in your event logger file. 
Once this is done, you also should add a block of code in the ``event_logger_struct_dump``
and the ``event_logger_struct_restore`` such that we can also restart with 
the log files. And lastly the ``event_logger_your_event_log_event`` function should
be called when your event takes place. If this is done than the only thing 
that needs to be done is including your struct in the ``event_logger_struct.h`` and 
you should be able to run with your new event logging feature. 

Adding a event logging debug mode
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Besides a global logging system in SWIFT that produces one file for the whole 
simulation, SWIFT also has the option to include debugging files for logging.
The debugging files are simpler to implement and logs every event
at the time it takes place. This means that the files become bigger because it 
contains all the individual events and it also produces a single file for every
node it runs on, so running on 4 nodes gives 4 files. There are only three
debugging functions: ``event_logger_your_event_init_log_file_debug`` which 
initializes your logger file; ``event_logger_your_event_init_debug`` which 
initializes your debugging struct (which contains only a lock and a file pointer);
``event_logger_your_event_log_event_debug`` which logs your event and locks the
file pointer. 

In the case of debugging you have access to more information than when you 
use the regular logger, this includes for example the individual particle 
ID(s) and the exact time step and times of the event. This allowes you to 
check if both are consistent with each other to verify that the global logging
system is working as implemented.

Implemented event logging
~~~~~~~~~~~~~~~~~~~~~~~~~

Currently there are 3 types of event logging implemented:

1. SNII logging (including a debugging mode)

2. SNII logging (including a debugging mode)

3. R-processes logging (no debugging mode)

These implemented event logging modes need to be specified in the parameter file
like:

.. code:: YAML

    # Event logging properties 
    Event_logger:
      delta_time_SNIa_Myr:        10
      delta_time_SNII_Myr:        1
      delta_time_r_processes_Myr: 20

