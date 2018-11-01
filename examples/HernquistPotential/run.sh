#!/bin/bash

# Generate the initial conditions if they are not present.
if [ ! -e Hernquist.hdf5 ]
then
    echo "Generating initial conditions for the isothermal potential box example..."
    python3 makeIC.py 1000 0
fi

rm -rf hernquist_*.hdf5
../swift -g -t 1 hernquist.yml 2>&1 | tee output.log


