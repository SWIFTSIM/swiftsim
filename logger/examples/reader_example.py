#!/usr/bin/env python3
"""
Read a logger file by using an index file.
Example: ./reader_example.py ../../examples/SedovBlast_3D/index 0.1
"""
import sys
import numpy as np
import matplotlib.pyplot as plt
from glob import glob
from mpl_toolkits.mplot3d import Axes3D
sys.path.append("../.libs/")

import liblogger as logger


def plot3D(pos):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    ax.plot(pos[:, 0], pos[:, 1], pos[:, 2], ".", markersize=0.1)


def plot2D(pos, entropies):
    plt.figure()
    center = np.array([0.5]*3)
    r2 = np.sum((pos - center)**2, axis=1)

    # plot entropy vs distance
    plt.plot(np.sqrt(r2), entropies, '.',
             markersize=0.2)

    plt.xlabel("Radius")
    plt.ylabel("Entropy")


filenames = "../../examples/HydroTests/SedovBlast_3D/index_*dump"
filenames = glob(filenames)
time = 0.01
if len(sys.argv) >= 2:
    time = float(sys.argv[1])
else:
    print("No time supplied (second argument), using default.")
if len(sys.argv) >= 3:
    filenames = sys.argv[2:]
else:
    print("No filenames supplied (first argument), using default.")

print("basename: %s" % filenames)
print("time: %g" % time)

# read the logger
pos = None
ent = None
for f in filenames:
    filename = f[:-5]
    with logger.Reader(filename, verbose=0) as reader:
        t = reader.get_time_limits()
        pos_tmp, ent_tmp = reader.get_particle_data(
            ["Coordinates", "Entropies"], time)

        # add the data to the list
        if pos is None:
            pos = pos_tmp
            ent = ent_tmp
        else:
            pos = np.append(pos, pos_tmp, axis=0)
            ent = np.append(ent, ent_tmp, axis=0)

print("Min/Max of the position:", pos.min(), pos.max())
print("Min/Max of the entropy:", ent.min(), ent.max())
plot3D(pos)
plot2D(pos, ent)
plt.show()
