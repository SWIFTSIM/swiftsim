###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2019 Loic Hausammann (loic.hausammann@epfl.ch)
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

import h5py
import numpy as np
from scipy.optimize import bisect

# Generates a swift IC file for the Rayleigh-Taylor vortex in a periodic box
# following the conditions given in Hopkins 2013

# Parameters
N = [50, 100]  # Particles along one edge
gamma = 5./3.  # Gas adiabatic index
delta = 0.025  # interface size
dv = 0.025  # velocity perturbation
rho1 = 1    # Lower region density
rho2 = 2    # Upper region density
P0 = 0      # integration constant of the hydrostatic equation
g = -0.5    # gravitational acceleration
box_size = [0.5, 1.]  # size of the box
fileOutputName = "rayleigh_taylor.hdf5"


def density(y):
    """
    Mass density as function of the position y.
    """
    tmp = 1 + np.exp(-(y - 0.5 * box_size[1]) / delta)
    return rho1 + (rho2 - rho1) / tmp


def mass(y):
    """
    Integrated Mass in x and y.
    """
    tmp = np.log(np.exp(y / delta) + np.exp(box_size[1] * 0.5 / delta))
    tmp -= 0.5 * box_size[1] / delta
    m = box_size[0] * (rho1 * y + (rho2 - rho1) * delta * tmp)
    return m


def inv_mass(m):
    """
    Inverse of the function `mass`.
    """
    left = 0
    right = box_size[1]

    def f(y, x):
        return x - mass(y)

    return bisect(f, left, right, args=(m))


def pressure(y):
    """
    Pressure as function of the position y.
    Here we assume hydrostatic equation.
    """
    tmp = np.log(np.exp(y / delta) + np.exp(box_size[1] * 0.5 / delta))
    tmp -= 0.5 * box_size[1] / delta
    tmp *= g * (rho2 - rho1) * delta
    return P0 - rho1 * g * y + tmp


# ---------------------------------------------------

if (box_size[0] / N[0] != box_size[1] / N[1]):
    raise Exception(
        "Assuming the same number of particle per unit of distance")


# Start by generating grids of particles
numPart = N[0] * N[1]

m_tot = mass(box_size[1])
m_part = m_tot / numPart

coords = np.zeros((numPart, 3))
h = np.ones(numPart) * 1.2348 * box_size[0] / N[0]
m = np.ones(numPart) * m_part
u = np.zeros(numPart)
vel = np.zeros((numPart, 3))


# generate grid of particles
y_prev = 0
# loop over y
for j in range(N[1]):
    m_y = m_part * N[0] + mass(y_prev)
    if m_y > m_tot:
        y_j = box_size[1] - 1e-5 * (box_size[1] - y_prev)
    else:
        y_j = inv_mass(m_y)
    y_prev = y_j

    # loop over x
    for i in range(N[0]):

        index = j * N[0] + i

        x = i * box_size[0] / float(N[0])

        coords[index, 0] = x
        coords[index, 1] = y_j
        u[index] = pressure(y_j) / (rho1 * (gamma-1.))

ids = np.linspace(1, numPart, numPart)

# Velocity perturbation
ind = coords[:, 1] < 0.7
ind = np.logical_and(ind, coords[:, 1] > 0.3)
vel[ind, 1] = (1 + np.cos(8 * np.pi * (coords[index, 0] + 0.5 * box_size[0])))
vel[ind, 1] += (1 + np.cos(5 * np.pi * (coords[index, 1] - 0.5 * box_size[1])))
vel[ind, 1] *= dv

# File
fileOutput = h5py.File(fileOutputName, 'w')

# Header
grp = fileOutput.create_group("/Header")
grp.attrs["BoxSize"] = box_size
grp.attrs["NumPart_Total"] = [numPart, 0, 0, 0, 0, 0]
grp.attrs["NumPart_Total_HighWord"] = [0, 0, 0, 0, 0, 0]
grp.attrs["NumPart_ThisFile"] = [numPart, 0, 0, 0, 0, 0]
grp.attrs["Time"] = 0.0
grp.attrs["NumFileOutputsPerSnapshot"] = 1
grp.attrs["MassTable"] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
grp.attrs["Flag_Entropy_ICs"] = [0, 0, 0, 0, 0, 0]
grp.attrs["Dimension"] = 2

# Units
grp = fileOutput.create_group("/Units")
grp.attrs["Unit length in cgs (U_L)"] = 1.
grp.attrs["Unit mass in cgs (U_M)"] = 1.
grp.attrs["Unit time in cgs (U_t)"] = 1.
grp.attrs["Unit current in cgs (U_I)"] = 1.
grp.attrs["Unit temperature in cgs (U_T)"] = 1.

# Particle group
grp = fileOutput.create_group("/PartType0")
ds = grp.create_dataset('Coordinates', (numPart, 3), 'd')
ds[()] = coords
ds = grp.create_dataset('Velocities', (numPart, 3), 'f')
ds[()] = vel
ds = grp.create_dataset('Masses', (numPart, 1), 'f')
ds[()] = m.reshape((numPart, 1))
ds = grp.create_dataset('SmoothingLength', (numPart, 1), 'f')
ds[()] = h.reshape((numPart, 1))
ds = grp.create_dataset('InternalEnergy', (numPart, 1), 'f')
ds[()] = u.reshape((numPart, 1))
ds = grp.create_dataset('ParticleIDs', (numPart, 1), 'L')
ds[()] = ids.reshape((numPart, 1))

fileOutput.close()
