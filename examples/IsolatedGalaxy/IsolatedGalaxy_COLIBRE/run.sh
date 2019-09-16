#!/bin/bash

if [ ! -e lowres64.hdf5 ] 
then     
    echo "Fetching initial conditions for the isolated galaxy example..."
    ./getIC.sh
fi

if [ ! -e yieldtables.tar.gz ] 
then     
    echo "Fetching yield tables for the isolated galaxy example..."
    ./getEagleYieldTable.sh
fi

../../swift --threads=8 --external-gravity --self-gravity --stars --star-formation --cooling --temperature --hydro --feedback --limiter isolated_galaxy.yml 2>&1 | tee output.log
