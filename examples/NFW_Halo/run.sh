#!/bin/bash

if [ ! -i test_nfw.hdf5 ]
    if command -v python3 &>/dev/null; then
        python3 makeIC.py
    else 
        python makeIC.py
    fi
fi

# self gravity G, external potential g, hydro s, threads t and high verbosity v
../swift -g -t 6 test.yml 2>&1 | tee output.log



if command -v python3 &>/dev/null; then
    python3 makePlots.py
else 
    python makePlots.py
fi
