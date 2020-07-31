#!/bin/bash

# filename=fid.hdf5
# filename=f10.hdf5
# filename=f90.hdf5
filename=lowres8.hdf5
# filename=lowres64.hdf5
# filename=lowres512.hdf5
# filename=highres6.hdf5

if [ ! -e $filename ]
then
    echo "Fetching initial conditons for the isolated galaxy with an external potential ..."
    ./getIC.sh $filename
fi

echo "Generate initial conditions"
python3 makeIC.py $filename

../../swift -v 2 --hydro --sinks --external-gravity --self-gravity --threads=16 isolated_galaxy.yml 2>&1 | tee output.log
