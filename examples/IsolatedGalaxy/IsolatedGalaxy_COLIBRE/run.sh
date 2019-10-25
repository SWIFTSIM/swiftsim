#!/bin/bash

if [ ! -e lowres64.hdf5 ] 
then     
    echo "Fetching initial conditions for the isolated galaxy example..."
    ./getIC.sh
fi

if [ ! -e UV_dust1_CR1_G1_shield1.hdf5 ]
then
    echo "Fetching cooling tables for the isolated galaxy example..."
    ./getColibreCoolingTables.sh
fi

if [ ! -e HIIregions_BPASS_binary.hdf5 ]
then
    echo "Fetching feedback tables for the isolated galaxy example..."
    ./getColibreFeedbackTables.sh
fi

if [ ! -e yieldtables.tar.gz ] 
then     
    echo "Fetching yield tables for the isolated galaxy example..."
    ./getEagleYieldTable.sh
fi

../../swift --threads=8 --external-gravity --self-gravity --stars --star-formation --cooling --temperature --hydro --feedback --limiter isolated_galaxy.yml 2>&1 | tee output.log
