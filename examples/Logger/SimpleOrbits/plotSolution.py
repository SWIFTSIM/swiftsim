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
id_focus = 0

G = 1.189972e-04
M = 3.329460e+05


def getSnapshots():
    filenames = glob("simple_orbits_*.hdf5")
    N = len(filenames)
    filenames.sort()

    E = np.zeros((N, makeIC.num_part))
    t = np.zeros(N)
    p = np.zeros((N, 3))

    for i, f in enumerate(filenames):
        if i % 10 == 0:
            print("File %i / %i" % (i, N))
        f = File(f, "r")
        pos = f["PartType1/Coordinates"][:]
        pos -= center
        vel = f["PartType1/Velocities"][:]
        ids = f["PartType1/ParticleIDs"][:]

        t[i] = f["Header"].attrs["Time"]
        for j in range(makeIC.num_part):
            ind = ids == j
            r = np.sum(pos[ind, :]**2)**0.5
            v2 = np.sum(vel[ind, :]**2)
            E[i, ind] = 0.5 * v2 - G * M / r
            if j == id_focus:
                p[i, :] = pos[ind, :]

    # compute the plotting variables
    eps = (E - E[0, :]) / E[0, :]
    return t, eps, p


def getStatistics():
    data = np.genfromtxt("energy.txt", names=True)

    times = data["Time"]
    E = data["E_tot"]
    eps = (E - E[0]) / E[0]
    return times, eps


def getLogger():
    basename = "index"
    N = 1000
    verbose = 0

    # Get time limits
    t_min, t_max = logger.getTimeLimits(basename, verbose)
    times = np.linspace(t_min, t_max, N)

    # Create output arrays
    E = np.zeros((N, makeIC.num_part))
    E_parts = np.zeros((N, makeIC.num_part))
    p = np.zeros((N, 3))
    t_parts = np.zeros((N, makeIC.num_part))

    # Read the particles
    parts = logger.loadSnapshotAtTime(
        basename, times[0], verbose)

    for i, t in enumerate(times):
        interp = logger.moveForwardInTime(
            basename, parts, t, verbose)
        rel_pos = interp["positions"] - center
        v = interp["velocities"]

        rel_pos_parts = parts["positions"] - center
        v_parts = parts["velocities"]

        for j in range(makeIC.num_part):
            # do interp
            ind = interp["ids"] == j
            r = np.sum(rel_pos[ind, :]**2)**0.5
            v2 = np.sum(v[ind, :]**2)
            E[i, ind] = 0.5 * v2 - G * M / r
            if j == id_focus:
                p[i, :] = rel_pos[ind, :]

            # do parts
            ind = parts["ids"] == j
            r = np.sum(rel_pos_parts[ind, :]**2)**0.5
            v2 = np.sum(v_parts[ind, :]**2)
            E_parts[i, ind] = 0.5 * v2 - G * M / r
            t_parts[i, ind] = parts["times"][ind]

    # compute the plotting variables
    eps = (E - E[0, :]) / E[0, :]
    eps_parts = (E_parts - E_parts[0, :]) / E_parts[0, :]
    return times, eps, p, eps_parts, t_parts


t_log, eps_log, p_log, eps_parts, t_parts = getLogger()
t_snap, eps_snap, p_snap = getSnapshots()
t_stat, eps_stat = getStatistics()

plt.figure()
plt.plot(t_snap, eps_snap, ".", label="Snapshot")
plt.plot(t_parts, eps_parts, "x", label="Logger")
plt.plot(t_log, eps_log, "--", label="Logger (Interpolation)")
plt.plot(t_stat, eps_stat, "-", label="Statistics")
plt.xlabel("Time [yr]")
plt.ylabel("Relative error on the specific energy")
plt.legend(ncol=2)


plt.figure()
plt.plot(p_snap[:, 0], p_snap[:, 1], ".", label="Snapshot")
plt.plot(p_log[:, 0], p_log[:, 1], "--", label="Logger")
plt.xlabel("Position [AU]")
plt.ylabel("Position [AU]")
plt.legend()
plt.show()
