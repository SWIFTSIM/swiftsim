#!/bin/bash

if [ ! -e lowres64.hdf5 ] 
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

../../swift --threads=16 --feedback --external-gravity --self-gravity --stars --star-formation --cooling --hydro isolated_galaxy.yml -P SNIaDTD:SNIa_efficiency_p_Msun:0.002 2>&1 | tee output.log

python3 plotSolution.py nu0.002

../../swift --threads=16 --feedback --external-gravity --self-gravity --stars --star-formation --cooling --hydro isolated_galaxy.yml -P SNIaDTD:SNIa_efficiency_p_Msun:0.001 2>&1 | tee output.log

python3 plotSolution.py nu0.001

../../swift --threads=16 --feedback --external-gravity --self-gravity --stars --star-formation --cooling --hydro isolated_galaxy.yml -P SNIaDTD:SNIa_efficiency_p_Msun:0.004 2>&1 | tee output.log

python3 plotSolution.py nu0.004

../../swift --threads=16 --feedback --external-gravity --self-gravity --stars --star-formation --cooling --hydro isolated_galaxy.yml -P SNIaDTD:normalization_timescale_Gyr:6.0 2>&1 | tee output.log

python3 plotSolution.py nu0.002_norm6

../../swift --threads=16 --feedback --external-gravity --self-gravity --stars --star-formation --cooling --hydro isolated_galaxy.yml -P SNIaDTD:SNIa_delay_time_Gyr:0.08 2>&1 | tee output.log

python3 plotSolution.py nu0.002_t_init
