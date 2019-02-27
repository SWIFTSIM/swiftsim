#!/bin/bash

beta=(1.0 0.1 0.01 0.001)
swift_location="../../swift"
flags="--hydro --limiter --threads=2"
parameter="diffusion.yml"
make_ic_script="makeIC.py"
plot_script="plotSolution.py"

for i in ${beta[@]};
do
    mkdir beta_$i
    cd beta_$i

    python ../$make_ic_script

    ../$swift_location $flags -P "SPH:diffusion_beta:${i}" ../$parameter

    python ../$plot_script

    rm *.hdf5

    cd ..
done
