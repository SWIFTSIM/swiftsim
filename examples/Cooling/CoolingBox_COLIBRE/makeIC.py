###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2016 Stefan Arridge (stefan.arridge@durhama.ac.uk)
#                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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

# Generates a SWIFT IC file with a constant density and pressure

# Parameters
periodic = 1       # 1 For periodic box
boxSize = 1.0      # kpc 
nH = 0.01          # Hydrogen atoms per cm^3 
T = 1.0e7          # Initial Temperature
gamma = 5./3.      # Gas adiabatic index
XH = 0.73738788833 # Hydrogen mass fraction. 
fileName = "coolingBox.hdf5"
# ---------------------------------------------------

# Code units in cgs. These must agree 
# with the units defined in the 
# parameter file. 
UnitMass_in_cgs = 1.9891E43        # 10^10 solar masses 
UnitLength_in_cgs = 3.08567758E21  # 1 kpc
UnitVelocity_in_cgs = 1.0e5        # km/s 

# Derived units. 
UnitDensity_in_cgs = UnitMass_in_cgs / (UnitLength_in_cgs ** 3.0) 
UnitTime_in_cgs = UnitLength_in_cgs / UnitVelocity_in_cgs 

# Physical constants 
m_h_cgs = 1.67e-24
k_b_cgs = 1.38e-16

# defines some constants
# need to be changed in plotTemperature.py too
h_frac = 0.76

# Read id, position and h from glass
glass = h5py.File("glassCube_32.hdf5", "r")
ids = glass["/PartType0/ParticleIDs"][:]
pos = glass["/PartType0/Coordinates"][:, :] * boxSize
h = glass["/PartType0/SmoothingLength"][:] * boxSize

# Compute basic properties
numPart = np.size(pos) // 3

# Density in cgs units 
rho_cgs = m_h_cgs * nH / XH 

# Density in code units 
rho = rho_cgs / UnitDensity_in_cgs 

# Particle mass in code units 
mass = boxSize**3 * rho / numPart

# Use the mean molecular weight for 
# either ionised or neutral gas, 
# based on the temperature. 
mu_neutral = 4. / (1. + 3. * XH)
mu_ionised = 4. / (8. - 5. * (1. - XH));
if (T > 1.0e4): 
    mu = mu_ionised 
else: 
    mu = mu_neutral 

internalEnergy = k_b_cgs * T * mu / ((gamma - 1.) * m_h_cgs)
internalEnergy /= UnitVelocity_in_cgs ** 2.0 

# File
f = h5py.File(fileName, 'w')

# Header
grp = f.create_group("/Header")
grp.attrs["BoxSize"] = boxSize
grp.attrs["NumPart_Total"] = [numPart, 0, 0, 0, 0, 0]
grp.attrs["NumPart_Total_HighWord"] = [0, 0, 0, 0, 0, 0]
grp.attrs["NumPart_ThisFile"] = [numPart, 0, 0, 0, 0, 0]
grp.attrs["Time"] = 0.0
grp.attrs["NumFilesPerSnapshot"] = 1
grp.attrs["MassTable"] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
grp.attrs["Flag_Entropy_ICs"] = 0

# Runtime parameters
grp = f.create_group("/RuntimePars")
grp.attrs["PeriodicBoundariesOn"] = periodic

# Units
grp = f.create_group("/Units")
grp.attrs["Unit length in cgs (U_L)"] = UnitLength_in_cgs
grp.attrs["Unit mass in cgs (U_M)"] = UnitMass_in_cgs
grp.attrs["Unit time in cgs (U_t)"] = UnitTime_in_cgs
grp.attrs["Unit current in cgs (U_I)"] = 1.
grp.attrs["Unit temperature in cgs (U_T)"] = 1.

# Particle group
grp = f.create_group("/PartType0")

v = np.zeros((numPart, 3))
ds = grp.create_dataset('Velocities', (numPart, 3), 'f')
ds[()] = v

m = np.full((numPart, 1), mass)
ds = grp.create_dataset('Masses', (numPart, 1), 'f')
ds[()] = m

h = np.reshape(h, (numPart, 1))
ds = grp.create_dataset('SmoothingLength', (numPart, 1), 'f')
ds[()] = h

u = np.full((numPart, 1), internalEnergy)
ds = grp.create_dataset('InternalEnergy', (numPart, 1), 'f')
ds[()] = u

ids = np.reshape(ids, (numPart, 1))
ds = grp.create_dataset('ParticleIDs', (numPart, 1), 'L')
ds[()] = ids

ds = grp.create_dataset('Coordinates', (numPart, 3), 'd')
ds[()] = pos

f.close()

print("Initial condition generated")
