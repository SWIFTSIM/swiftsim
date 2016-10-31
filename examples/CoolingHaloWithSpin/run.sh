#!/bin/bash

# Generate the initial conditions if they are not present.
echo "Generating initial conditions for the isothermal potential box example..."
python makeIC.py 1000000

../swift -g -s -C -D -t 16 cooling_halo.yml 2>&1 | tee output.log

python radial_profile.py 2. 100 101

python internal_energy_profile.py 2. 100 101

python test_energy_conservation.py
