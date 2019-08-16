#!/bin/bash

if [ ! -e lowres64.hdf5 ]
then
echo "Fetching initial conditions for the isolated galaxy example..."
./getIC.sh
fi

../../swift --threads=16 --feedback --external-gravity --self-gravity --stars --star-formation --cooling --hydro isolated_galaxy.yml 2>&1 | tee output.log

