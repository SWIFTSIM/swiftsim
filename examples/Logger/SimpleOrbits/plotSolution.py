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
from scipy.integrate import solve_ivp
import makeIC
import sys
sys.path.append("../../../logger/.libs/")
import liblogger as logger

# Plot parameters
plt.style.use("mpl_style")

center = np.array([0.5 * makeIC.boxsize]*3)
id_focus = 0
# Plot the absolute positions (+vel) or the difference
# to a solution computed with a RK solver.
diff_to_sol = False

# Defines the constants
G = 1.189972e-04  # gravitational constant
M = 3.329460e+05  # mass of the sun

# generate the figures
fig_1 = plt.figure()
fig_2 = plt.figure()
fig_3 = plt.figure()


def gravity(t, y):
    """
    Compute the equation of motion.

    Parameters
    ----------

    t: float
      Time of the step
    y: np.array
      Variable to evolve [x, y, vx, vy].
    """
    dy = np.zeros(4)
    dy[:2] = y[2:]

    r = np.sum((y[:2] - center[:2])**2)**0.5

    a = - G * M / r**3
    dy[2] = y[2] * a
    dy[3] = y[3] * a
    return dy


def Evolve(y0, t):
    """
    Compute the solution to the gravitational problem.

    Parameters
    ----------

    y0: np.array
      Initial conditions [x, y, vx, vy]

    t: np.array
      Time of the outputs.
    """
    y = solve_ivp(gravity, (t[0], t[-1]), y0, t_eval=t)
    return y.y


def plotRelative(t, p, *args, **kwargs):
    """
    Wrapper around the function plot from matplotlib, but
    plot the relative evolution of the variable.
    """
    p = (p - p[0]) / p[0]
    plt.plot(t, p, *args, **kwargs)


def doSnapshots():
    """
    Read the snapshots and plot the corresponding variables.
    """

    # get all the filenames
    filenames = glob("simple_orbits_*.hdf5")
    N = len(filenames)
    filenames.sort()

    # generate the output arrays
    E = np.zeros((N, makeIC.num_part))
    t = np.zeros(N)
    p = np.zeros((N, 3))
    v = np.zeros((N, 3))

    for i, f in enumerate(filenames):
        # get the data from the file
        f = File(f, "r")
        pos = f["PartType1/Coordinates"][:]
        pos -= center
        vel = f["PartType1/Velocities"][:]
        ids = f["PartType1/ParticleIDs"][:]

        t[i] = f["Header"].attrs["Time"]

        # loop over all the particles
        # (avoids comparing different particles)
        for j in range(makeIC.num_part):
            # Get the current id
            ind = ids == j

            # Compute the energy
            r = np.sum(pos[ind, :]**2)**0.5
            v2 = np.sum(vel[ind, :]**2)
            E[i, ind] = 0.5 * v2 - G * M / r

            # Get the pos / vel of the required particle
            if j == id_focus:
                p[i, :] = pos[ind, :]
                v[i, :] = vel[ind, :]

    # Compute the solution
    y0 = np.zeros(4)
    y0[:2] = p[0, :2]
    y0[2:] = v[0, :2]
    y = Evolve(y0, t)

    # compute the plotting variables
    plt.figure(fig_1.number)
    plotRelative(t, E, ".", label="Snapshot")

    plt.figure(fig_2.number)
    if diff_to_sol:
        p[:, :2] = p[:, :2] - y[:2, :].transpose()
    plt.plot(p[:, 0], p[:, 1], ".", label="Snapshot")

    plt.figure(fig_3.number)
    if diff_to_sol:
        v[:, :2] = v[:, :2] - y[2:, :].transpose()
    plt.plot(v[:, 0], v[:, 1], ".", label="Snapshot")


def doStatistics():
    """
    Do the plots with the energy output file.
    """
    data = np.genfromtxt("energy.txt", names=True)

    times = data["Time"]
    E = data["E_tot"]
    plt.figure(fig_1.number)
    plotRelative(times, E, "-", label="Statistics")


def doLogger():
    """
    Read the logfile and plot the corresponding variables.
    """
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
    v = np.zeros((N, 3))
    t_parts = np.zeros((N, makeIC.num_part))
    p_parts = np.zeros((N, 3))
    v_parts = np.zeros((N, 3))

    # Read the particles
    parts = logger.loadSnapshotAtTime(
        basename, times[0], verbose)

    for i, t in enumerate(times):
        # Get the next particles
        interp = logger.moveForwardInTime(
            basename, parts, t, verbose)
        rel_pos = interp["positions"] - center
        vel = interp["velocities"]

        rel_pos_parts = parts["positions"] - center
        vel_parts = parts["velocities"]

        for j in range(makeIC.num_part):
            # Compute the interpolated variables
            ind = interp["ids"] == j
            r = np.sum(rel_pos[ind, :]**2)**0.5
            v2 = np.sum(vel[ind, :]**2)
            E[i, ind] = 0.5 * v2 - G * M / r
            if j == id_focus:
                p[i, :] = rel_pos[ind, :]
                v[i, :] = vel[ind, :]

            # Compute the variables of the last record
            ind = parts["ids"] == j
            r = np.sum(rel_pos_parts[ind, :]**2)**0.5
            v2 = np.sum(vel_parts[ind, :]**2)
            E_parts[i, ind] = 0.5 * v2 - G * M / r
            t_parts[i, ind] = parts["times"][ind]
            if j == id_focus:
                p_parts[i, :] = rel_pos_parts[ind, :]
                v_parts[i, :] = vel_parts[ind, :]

    # compute the plotting variables
    plt.figure(fig_1.number)
    plotRelative(t_parts, E_parts, "x", label="Logger")
    plotRelative(times, E, "--", label="Logger (Interpolation)")

    # Compute the solution
    y0 = np.zeros(4)
    y0[:2] = p[0, :2]
    y0[2:] = v[0, :2]
    t_parts, ind = np.unique(t_parts[:, 0], return_index=True)
    y = Evolve(y0, t_parts)
    if diff_to_sol:
        p_parts = p_parts[ind, :2] - y[:2, :].transpose()

    # plot the solution
    plt.figure(fig_2.number)
    plt.plot(p_parts[:, 0], p_parts[:, 1], "x", label="Logger")

    plt.figure(fig_3.number)
    if diff_to_sol:
        v_parts = v_parts[ind, :2] - y[2:, :].transpose()
    plt.plot(v_parts[:, 0], v_parts[:, 1], "x", label="Logger")

    # Compute the solution
    y0 = np.zeros(4)
    y0[:2] = p[0, :2]
    y0[2:] = v[0, :2]
    y = Evolve(y0, times)
    if diff_to_sol:
        p[:, :2] = p[:, :2] - y[:2, :].transpose()
    plt.figure(fig_2.number)
    plt.plot(p[:, 0], p[:, 1], "--", label="Logger (Interpolation)")

    plt.figure(fig_3.number)
    if diff_to_sol:
        v[:, :2] = v[:, :2] - y[2:, :].transpose()
    plt.plot(v[:, 0], v[:, 1], "--", label="Logger (Interpolation)")


# do all the plots
doStatistics()
doSnapshots()
doLogger()

# add text
plt.figure(fig_1.number)
plt.xlabel("Time [yr]")
plt.ylabel(r"$\frac{E - E(t=0)}{E(t=0)}$")
plt.legend(ncol=2)

plt.figure(fig_2.number)
plt.xlabel("Position [AU]")
plt.ylabel("Position [AU]")
plt.legend()

plt.figure(fig_3.number)
plt.xlabel("Velocity [AU / yr]")
plt.ylabel("Velocity [AU / yr]")
plt.legend()

plt.show()
