#!/bin/bash

if [ ! -e lowres64.hdf5 ] 
then     
    echo "Fetching initial conditions for the isolated galaxy example..."
    ./getIC64.sh
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

../../swift --threads=16 --feedback --external-gravity --self-gravity --stars --star-formation --cooling --hydro ../isolated_galaxy.yml -P EAGLEEntropyFloor:Jeans_temperature_norm_K:100. -P InitialConditions:file_name:lowres64.hdf5 -P Gravity:max_physical_baryon_softening:0.8 2>&1 | tee output.log

