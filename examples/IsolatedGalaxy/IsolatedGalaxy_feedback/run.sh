#!/bin/bash

if [ ! -e lowres8.hdf5 ] 
then     
    echo "Fetching initial conditions for the isolated galaxy example..."
    ./getIC.sh
fi

if [ ! -e coolingtables ] 
then     
    echo "Fetching EAGLE cooling tables for the isolated galaxy example..."
    ./getEagleCoolingTable.sh
fi

if [ ! -e yieldtables ] 
then     
    echo "Fetching EAGLE stellar yield tables for the isolated galaxy example..."
    ./getYieldTable.sh
fi

../../swift --threads=16 --feedback --external-gravity --self-gravity --stars --star-formation --cooling --hydro ../isolated_galaxy.yml -P InitialConditions:file_name:lowres8.hdf5 -P Gravity:max_physical_baryon_softening:0.4 2>&1 | tee output.log

# Kennicutt-Schmidt law plot
python3 ../plotKSlaw.py

# Plot that the random star formation matches the expected SFH
python3 ../SFH.py

# Plot of the box enrichment
python3 ../plot_box_evolution.py
