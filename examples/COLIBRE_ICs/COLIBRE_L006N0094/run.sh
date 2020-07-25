#!/bin/bash

# Generate the initial conditions if they are not present.
if [ ! -e Eagle_06Mpc_094.hdf5 ]
then
    echo "Fetching initial conditions for the COLIBRE 12Mpc example..."
    ./getIC.sh
fi

# Grab the cooling and yield tables if they are not present.
if [ ! -e yieldtables ]
then
    echo "Fetching EAGLE yield tables..."
    ../getEagleYieldTable.sh
fi

if [ ! -e UV_dust1_CR1_G1_shield1.hdf5 ]
then
    echo "Fetching COLIBRE cooling tables..."
    ../getColibreCoolingTables.sh
fi

if [ ! -e HIIregions_BPASS_binary.hdf5 ]
then
    echo "Fetching COLIBRE early feedback tables..."
    ../getColibreFeedbackTables.sh
fi

# The following run-time options are broken down by line as:
# Basic run-time options
# Create and run with stars
# Radiative options - run with cooling and stellar feedback
# Run with the time-step limiter required to capture feedback
# Run with black holes - fof is needed for the seeding
# Threading options - run with threads and pinning (latter not required but improves performance)
# The corresponding parameter file for this run

../../swift \
    --cosmology --eagle \
    --threads=16 --pin \
    colibre_06Mpc_094.yml

