#!/bin/bash
rm -f snapshots/uranus_impact_??????.hdf5
rm -f timesteps_*.txt

../swift -G -s -t 8 uranus_impact.yml
