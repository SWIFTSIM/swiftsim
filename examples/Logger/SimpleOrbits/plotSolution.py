###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2020 Loic Hausammann (loic.hausammann@epfl.ch)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
##############################################################################

from h5py import File
import numpy as np
from glob import glob
import matplotlib.pyplot as plt
import makeIC
import sys
sys.path.append("../../../logger/.libs/")
import liblogger as logger

# Plot parameters
plt.style.use("mpl_style")

center = np.array([0.5 * makeIC.boxsize]*3)


def getErrorSnapshots():
    filenames = glob("simple_orbits_*.hdf5")
    N = len(filenames)
    filenames.sort()

    r2 = np.zeros((N, makeIC.num_part))
    t = np.zeros(N)

    for i, f in enumerate(filenames):
        if i % 10 == 0:
            print("File %i / %i" % (i, N))
        f = File(f, "r")
        pos = f["PartType1/Coordinates"][:]
        pos -= center
        ids = f["PartType1/ParticleIDs"][:]
        print(ids)

        t[i] = f["Header"].attrs["Time"]
        for j in range(makeIC.num_part):
            ind = ids == j
            print(ind)
            r2[i, j] = np.sum(pos[ind, :]**2)

    # compute the plotting variables
    r = np.sqrt(r2)
    av = r.mean(axis=0)
    eps = (r - av) / r[0, :]
    return t, eps


def getErrorLogger():
    basename = "index"
    N = 1000
    t_min, t_max = logger.getTimeLimits(basename)
    times = np.linspace(t_min, t_max, N)
    r2 = np.zeros((N, makeIC.num_part))

    for i, t in enumerate(times):
        parts = logger.loadSnapshotAtTime(basename, t)
        pos = parts["positions"]
        pos -= center

        for j in range(makeIC.num_part):
            ind = parts["ids"] == j
            r2[i, j] = np.sum(pos[ind, :]**2)

    # compute the plotting variables
    r = np.sqrt(r2)
    av = r.mean(axis=0)
    eps = (r - av) / r[0, :]
    return times, eps


t_log, eps_log = getErrorLogger()
t_snap, eps_snap = getErrorSnapshots()

plt.figure()
colors = ["b", "r", "m", "g", "k"]
for i in range(makeIC.num_part):
    plt.plot(t_snap, eps_snap[:, i], ":", color=colors[i])
    plt.plot(t_log, eps_log[:, i], "-", color=colors[i])
plt.xlabel("Time [yr]")
plt.ylabel("Relative error on the radius")
plt.legend(["Snapshot", "Logger"])
plt.show()
