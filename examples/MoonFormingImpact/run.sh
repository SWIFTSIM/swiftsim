#!/bin/bash
rm -f snapshots/moon_forming_impact_??????.hdf5
rm -f timesteps_*.txt

../swift -G -s -t 8 moon_forming_impact.yml
