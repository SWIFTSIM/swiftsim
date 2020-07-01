#!/usr/bin/env python3
"""
Read a logger file by using an index file.
Example: ./reader_example.py ../../examples/SedovBlast_3D/index 0.1
"""
import sys
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
sys.path.append("../.libs/")

import liblogger as logger


def plot3D(data):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    if "gas" in data:
        pos = data["gas"]["Coordinates"]
        ax.plot(pos[:, 0], pos[:, 1], pos[:, 2], ".", markersize=0.1)
    if "dark_matter" in data:
        pos = data["dark_matter"]["Coordinates"]
        ax.plot(pos[:, 0], pos[:, 1], pos[:, 2], ".", markersize=0.1,
                color="k")


def plot2D(data):
    pos = data["gas"]["Coordinates"]

    center = np.array([0.5]*3)
    r2 = np.sum((pos - center)**2, axis=1)

    # plot entropy vs distance
    plt.plot(np.sqrt(r2), data["gas"]["Entropies"], '.',
             markersize=0.2)

    plt.xlim(0., 0.5)
    plt.ylim(-1, 50)
    plt.xlabel("Radius")
    plt.ylabel("Entropy")


basename = "../../examples/HydroTests/SedovBlast_3D/index_0000"
time = 0.03
if len(sys.argv) >= 2:
    basename = sys.argv[1]
else:
    print("No basename supplied (first argument), using default.")
if len(sys.argv) >= 3:
    time = float(sys.argv[2])
else:
    print("No time supplied (second argument), using default.")
if len(sys.argv) > 3:
    print("Ignoring excess arguments '%s'." % sys.argv[3:])
print("basename: %s" % basename)
print("time: %g" % time)

# read the logger
t = logger.getTimeLimits(basename)
data = logger.loadSnapshotAtTime(basename, time)

if "gas" in data:
    print("The data contains the following elements for the gas:")
    print(data["gas"].dtype.names)

if "dark_matter" in data:
    print("The data contains the following elements for the dark matter:")
    print(data["dark_matter"].dtype.names)

if "stars" in data:
    print("The data contains the following elements for the stars:")
    print(data["stars"].dtype.names)

plot3D(data)
# plot2D(data)
plt.show()
