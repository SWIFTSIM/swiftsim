#!/bin/bash

# Generate the initial conditions if they are not present.
if [ ! -e test_tensors.hdf5 ]
then
    echo "Generating initial conditions..."
    python makeIC.py
fi

rm -rf output_*.hdf5
../../swift --self-gravity --stars --threads=1 test.yml 2>&1 | tee output.log

python print_tensors.py output_0000.hdf5
