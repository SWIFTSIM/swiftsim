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
    time = 0.05
else:
    basename = sys.argv[-1]
    time = sys.argv[-2]

# read dump
data = logger.loadFromIndex(basename, time)

pos = data["positions"]
print(pos.shape)


def plot3D():
    from mpl_toolkits.mplot3d import Axes3D
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    ax.plot(data["positions_x"], data["positions_y"],
            data["positions_z"], ".")


def plot2D():
    center = np.array([0.5]*3)
    r2 = np.sum((pos - center)**2, axis=1)

    # plot entropy vs distance
    plt.plot(np.sqrt(r2), data["entropies"], '.',
             markersize=0.2)

    plt.xlim(0., 0.5)
    plt.ylim(-1, 50)
    plt.xlabel("Radius")
    plt.ylabel("Entropy")


plot2D()
plt.show()
