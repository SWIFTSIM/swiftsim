#!/bin/bash

# MUSIC binary
music=~/music/Music
if test -f $music; then
    echo "Using the following version of MUSIC $music."
else
    echo "MUSIC is not found."
fi

echo "Generating the initial conditions"
$music music.conf

echo "Converting the initial conditions into a SWIFT compatible format"
python3 convert_ic.py

echo "Running SWIFT"
../../swift --cooling --feedback --cosmology  --limiter --sync --self-gravity --hydro --stars --star-formation --threads=24 zoom_in.yml 2>&1 | tee output.log
