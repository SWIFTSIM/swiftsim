#!/bin/bash

# Generate the initial conditions if they are not present.
if [ ! -e Hernquist.hdf5 ]
then
    echo "Generating initial conditions for the isothermal potential box example..."
    if command -v python3 &>/dev/null; then
        python3 makeIC.py 
    else 
        python makeIC.py
    fi
fi

rm -rf hernquist_*.hdf5
../swift -g -t 1 hernquist.yml 2>&1 | tee output.log



echo "Make plots of the radially free falling particles" 
if command -v python3 &>/dev/null; then
    python3 plotprog.py 
else 
    python plotprog.py
fi
