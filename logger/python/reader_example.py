#!/usr/bin/env python3
"""
Read a logger file by using an index.
Example: ./reader_example.py ../../examples/SedovBlast_3D/index 0.1
"""
import sys
import numpy as np
import matplotlib.pyplot as plt
sys.path.append("../.libs/")

import liblogger as logger

# Get filenames
if len(sys.argv) != 3:
    print("WARNING missing arguments. Will use the default ones")
    basename = "../../examples/HydroTests/SedovBlast_3D/index"
    time = 0.049
else:
    basename = sys.argv[-1]
    time = sys.argv[-2]

# read dump
data = logger.loadFromIndex(basename, time)

# Compute distance from center
plt.show()
pos = np.array([data["positions_x"], data["positions_y"],
               data["positions_y"]]).transpose()
center = pos.mean(axis=0)
r2 = np.sum((pos - center)**2, axis=1)

# plot entropy vs distance
plt.plot(np.sqrt(r2), data["densities"], '.')

plt.xlim(0., 0.5)
plt.ylim(-1, 5)
plt.xlabel("Radius")
plt.ylabel("Density")

plt.show()
