#!/bin/bash

if [ ! -e COLIBRE_fid.hdf5 ] 
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

../../swift --threads=32 --feedback --external-gravity --self-gravity --stars --star-formation --cooling --hydro isolated_galaxy.yml 2>&1 | tee output.log

python3 plotSolution.py 50
