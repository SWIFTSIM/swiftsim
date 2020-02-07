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
from astropy import units, constants

np.random.seed(50)

# parameters

filename = "simple_orbits.hdf5"
num_part = 3
masses = 1.
# If changed, need to update simple_orbits.yml
M = units.solMass.to("earthMass")
M = float(M)

min_r = 0.2
max_r = 5
boxsize = 2 * max_r

u_l = 1.49e13  # AU
u_m = 5.97e27  # Earth mass
u_v = 474047.  # AU / yr
u_t = u_l / u_v
G = 1.189972e-04  # grav. const.

# generate the coordinates
dist = np.random.rand(num_part) * (max_r - min_r) + min_r
angle = np.random.rand(num_part) * 2 * np.pi
# more easy to do it in 2D, therefore coords[:, 2] == 0
coords = np.zeros((num_part, 3))
coords[:, 0] = dist * np.sin(angle)
coords[:, 1] = dist * np.cos(angle)
coords += boxsize * 0.5

# generate the masses
m = np.ones(num_part) * masses

# generate the ids
ids = np.arange(num_part)

# generate the velocities
sign = np.random.rand(num_part)
sign[sign < 0.5] = -1
sign[sign >= 0.5] = 1

v = np.zeros((num_part, 3))
v[:, 0] = sign * np.sqrt(G * M / (dist * (1 + np.tan(angle)**2)))
v[:, 1] = - np.tan(angle) * v[:, 0]

# File
file = File(filename, 'w')

# Header
grp = file.create_group("/Header")
grp.attrs["BoxSize"] = [boxsize]*3
grp.attrs["NumPart_Total"] = [0, num_part, 0, 0, 0, 0]
grp.attrs["NumPart_Total_HighWord"] = [0, 0, 0, 0, 0, 0]
grp.attrs["NumPart_ThisFile"] = [0, num_part, 0, 0, 0, 0]
grp.attrs["Time"] = 0.0
grp.attrs["NumFilesPerSnapshot"] = 1
grp.attrs["MassTable"] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
grp.attrs["Flag_Entropy_ICs"] = 0
grp.attrs["Dimension"] = 3

# Units
grp = file.create_group("/Units")
grp.attrs["Unit length in cgs (U_L)"] = u_l
grp.attrs["Unit mass in cgs (U_M)"] = u_m
grp.attrs["Unit time in cgs (U_t)"] = u_t
grp.attrs["Unit current in cgs (U_I)"] = 1.
grp.attrs["Unit temperature in cgs (U_T)"] = 1.

# Particle group
grp = file.create_group("/PartType1")
grp.create_dataset('Coordinates', data=coords, dtype='d')
grp.create_dataset('Velocities', data=v, dtype='f')
grp.create_dataset('Masses', data=m, dtype='f')
grp.create_dataset('ParticleIDs', data=ids, dtype='L')

file.close()
