#!/bin/bash

# Generate the initial conditions if they are not present.
if [ ! -e glassCube_32.hdf5 ]
then
    echo "Fetching initial glass file for the Supernovae feedback example..."
    ./getGlass.sh
fi
if [ ! -e stellar_evolution.hdf5 ]
then
    echo "Generating initial conditions for the 3D stellar evolution example..."
    python makeIC.py
fi

# Get the Yield tables
if [ ! -e yieldtables ]
then
    echo "Fetching Yield tables..."
    ./getEagleYieldTable.sh
fi

# Get the solutions
if [ ! -e StellarEvolutionSolution ]
then
    echo "Fetching solutions ..."
    ./getSolutions.sh
fi

../../swift --external-gravity --feedback --stars --hydro --threads=4 stellar_winds.yml  2>&1 | tee output.log
