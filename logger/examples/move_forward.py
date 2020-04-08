#!/usr/bin/env python3
"""
Generates a small movie using the moveForwardInTime function.
Example: ./move_forward.py ../../examples/SedovBlast_3D/index 100
"""
import sys
import numpy as np
import matplotlib.pyplot as plt
import swiftsimio.visualisation as vis
from subprocess import call
import os
sys.path.append("../.libs/")

import liblogger as logger

resolution = 1024

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


def plot2D():
    center = np.array([0.5]*3)
    r2 = np.sum((pos - center)**2, axis=1)

    # plot entropy vs distance
    plt.plot(np.sqrt(r2), data["Entropies"], '.',
             markersize=0.2)

    plt.xlim(0., 0.5)
    plt.ylim(-1, 50)
    plt.xlabel("Radius")
    plt.ylabel("Entropy")


basename = "../../examples/HydroTests/SedovBlast_3D/index_0000"
N = 100
if len(sys.argv) >= 2:
    basename = sys.argv[1]
else:
    print("No basename supplied (first argument), using default.")
if len(sys.argv) >= 3:
    N = int(sys.argv[2])
else:
    print("No number of points supplied (second argument), using default.")
if len(sys.argv) > 2:
    print("Ignoring excess arguments '%s'." % sys.argv[2:])
print("basename: %s" % basename)
print("N: %i" % N)

# read the logger
t = logger.getTimeLimits(basename)
print("Time limits: [%g, %g]" % t)
data = logger.loadSnapshotAtTime(basename, t[0])

if "gas" in data:
    print("The data contains the following elements for the gas:")
    print(data["gas"].dtype.names)

if "dark_matter" in data:
    print("The data contains the following elements for the dark matter:")
    print(data["dark_matter"].dtype.names)

if "stars" in data:
    print("The data contains the following elements for the stars:")
    print(data["stars"].dtype.names)

pos = data["gas"]["Coordinates"]

# Now create the movie.
A = np.zeros(N)
times = np.linspace(t[0], t[1], N)
for i, time in enumerate(times):
    verbose = 0
    new_array = 1  # Do we wish a new array or interpolate in place
    interp = logger.moveForwardInTime(basename, data, time, verbose,
                                      new_array)
    # Get the arrays
    pos = interp["gas"]["Coordinates"]
    h = interp["gas"]["SmoothingLengths"]
    m = interp["gas"]["Masses"]

    # make the required transformations
    width = pos.max() - pos.min()
    pos /= width
    h /= width

    plt.figure()
    img = vis.scatter(pos[:, 0], pos[:, 1], m, h, resolution)
    img /= width**2
    plt.imshow(img)
    plt.savefig("img_%04i.png" % i)
    plt.close()

# Convert the images into a movie
convert = "ffmpeg -i img_%04d.png -y -vcodec libx264 "
convert += "-profile:v high444 -refs 16 -crf 0 "
convert += "-preset ultrafast movie.mp4"
call(convert, shell=True)


# cleanup
os.remove("img_*.png")
