#!/bin/bash

# Generate the initial conditions if they are not present.
echo "Generating initial conditions for the isothermal potential box example..."
python makeIC.py 100000 

../swift -g -G -s -C -t 16 cooling_halo.yml 2>&1 | tee output.log

# python radial_profile.py 10

# python internal_energy_profile.py 10

# python test_energy_conservation.py
