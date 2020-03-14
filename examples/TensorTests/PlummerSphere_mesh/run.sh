#!/bin/bash

# Generate the initial conditions if they are not present.
if [ ! -e plummer.dat ]
then
    echo "Generating plummer model..."
    source make_plummer.sh
fi
if [ ! -e plummer.hdf5 ]
then
    echo "Generating initial conditions..."
    python makeIC.py
fi

rm -rf output_*.hdf5
../../swift --self-gravity --stars --threads=1 params.yml 2>&1 | tee output.log

#python plot.py
