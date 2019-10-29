#!/bin/bash

# # Generate the initial conditions if they are not present.
# if [ ! -e glassCube_32.hdf5 ]
# then
#     echo "Fetching initial glass file for the Supernovae feedback example..."
#     ./getGlass.sh
# fi
# if [ ! -e stellar_evolution.hdf5 ]
# then
#     echo "Generating initial conditions for the 3D stellar evolution example..."
#     python makeIC.py
# fi

# # Get the Yield tables
# if [ ! -e yieldtables ]
# then
#     echo "Fetching Yield tables..."
#     ./getEagleYieldTable.sh
# fi

# # Get the solutions
# if [ ! -e StellarEvolutionSolution ]
# then
#     echo "Fetching solutions ..."
#     ./getSolutions.sh
# fi

../../swift  --temperature --feedback --stars --hydro --external-gravity --cooling --threads=8 stellar_evolution_dust.yml -P COLIBREChemistry:init_abundance_metal:0.08 -P COLIBREChemistry:init_abundance_Hydrogen:0.752 -P COLIBREChemistry:init_abundance_Helium:0.248 > output_nosputter.txt
../../swift  --temperature --feedback --stars --hydro --external-gravity --cooling --threads=8 stellar_evolution_dust_sputter.yml -P COLIBREChemistry:init_abundance_metal:0.08 -P COLIBREChemistry:init_abundance_Hydrogen:0.752 -P COLIBREChemistry:init_abundance_Helium:0.248 > output_sputter.txt

python dustevo.py

echo "Done Z=0.08"

../../swift  --temperature --feedback --stars --hydro --external-gravity --cooling --threads=8 stellar_evolution_dust.yml -P COLIBREChemistry:init_abundance_metal:0.04 -P COLIBREChemistry:init_abundance_Hydrogen:0.752 -P COLIBREChemistry:init_abundance_Helium:0.248 > output_nosputter.txt
../../swift  --temperature --feedback --stars --hydro --external-gravity --cooling --threads=8 stellar_evolution_dust_sputter.yml -P COLIBREChemistry:init_abundance_metal:0.04 -P COLIBREChemistry:init_abundance_Hydrogen:0.752 -P COLIBREChemistry:init_abundance_Helium:0.248 > output_sputter.txt

python dustevo.py

echo "Done Z=0.04"

../../swift  --temperature --feedback --stars --hydro --external-gravity --cooling --threads=8 stellar_evolution_dust.yml -P COLIBREChemistry:init_abundance_metal:0.01 -P COLIBREChemistry:init_abundance_Hydrogen:0.752 -P COLIBREChemistry:init_abundance_Helium:0.248 > output_nosputter.txt
../../swift  --temperature --feedback --stars --hydro --external-gravity --cooling --threads=8 stellar_evolution_dust_sputter.yml -P COLIBREChemistry:init_abundance_metal:0.01 -P COLIBREChemistry:init_abundance_Hydrogen:0.752 -P COLIBREChemistry:init_abundance_Helium:0.248 > output_sputter.txt

python dustevo.py

echo "Done Z=0.01"

../../swift  --temperature --feedback --stars --hydro --external-gravity --cooling --threads=8 stellar_evolution_dust.yml -P COLIBREChemistry:init_abundance_metal:0.001 -P COLIBREChemistry:init_abundance_Hydrogen:0.752 -P COLIBREChemistry:init_abundance_Helium:0.248 > output_nosputter.txt
../../swift  --temperature --feedback --stars --hydro --external-gravity --cooling --threads=8 stellar_evolution_dust_sputter.yml -P COLIBREChemistry:init_abundance_metal:0.001 -P COLIBREChemistry:init_abundance_Hydrogen:0.752 -P COLIBREChemistry:init_abundance_Helium:0.248 > output_sputter.txt

python dustevo.py

echo "Done Z=0.001"

../../swift  --temperature --feedback --stars --hydro --external-gravity --cooling --threads=8 stellar_evolution_dust.yml -P COLIBREChemistry:init_abundance_metal:0.0001 -P COLIBREChemistry:init_abundance_Hydrogen:0.752 -P COLIBREChemistry:init_abundance_Helium:0.248 > output_nosputter.txt
../../swift  --temperature --feedback --stars --hydro --external-gravity --cooling --threads=8 stellar_evolution_dust_sputter.yml -P COLIBREChemistry:init_abundance_metal:0.0001 -P COLIBREChemistry:init_abundance_Hydrogen:0.752 -P COLIBREChemistry:init_abundance_Helium:0.248 > output_sputter.txt

python dustevo.py

echo "Done Z=1e-4"
