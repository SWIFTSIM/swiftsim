###############################################################################
 # This file is part of SWIFT.
 # Copyright (c) 2016 John A. Regan (john.a.regan@durham.ac.uk)
 #                    Tom Theuns (tom.theuns@durham.ac.uk)
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
import sys
import numpy
import math
import random
import matplotlib.pyplot as plt

# Generates a disc-patch in hydrostatic equilibrium
# see Creasey, Theuns & Bower, 2013, for the equations:
# disc parameters are: surface density sigma
#                      scale height b
# density: rho(z) = (sigma/2b) sech^2(z/b)
# isothermal velocity dispersion = <v_z^2? = b pi G sigma
# grad potential  = 2 pi G sigma tanh(z/b)
# potential       = ln(cosh(z/b)) + const
# Dynamical time  = sqrt(b / (G sigma))
# to obtain the 1/ch^2(z/b) profile from a uniform profile (a glass, say, or a uniform random variable), note that, when integrating in z
# \int 0^z dz/ch^2(z) = tanh(z)-tanh(0) = \int_0^x dx = x (where the last integral refers to a uniform density distribution), so that z = atanh(x)

# physical constants in cgs
NEWTON_GRAVITY_CGS  = 6.672e-8
SOLAR_MASS_IN_CGS   = 1.9885e33
PARSEC_IN_CGS       = 3.0856776e18
PROTON_MASS_IN_CGS  = 1.6726231e24
YEAR_IN_CGS         = 3.154e+7

# choice of units
const_unit_length_in_cgs   =   (PARSEC_IN_CGS)
const_unit_mass_in_cgs     =   (SOLAR_MASS_IN_CGS)
const_unit_velocity_in_cgs =   (1e5)

print "UnitMass_in_cgs:     ", const_unit_mass_in_cgs 
print "UnitLength_in_cgs:   ", const_unit_length_in_cgs
print "UnitVelocity_in_cgs: ", const_unit_velocity_in_cgs


# parameters of potential
surface_density = 100. # surface density of all mass, which generates the gravitational potential
scale_height    = 100.
gamma           = 5./3.
fgas            = 0.1  # gas fraction

# derived units
const_unit_time_in_cgs = (const_unit_length_in_cgs / const_unit_velocity_in_cgs)
const_G                = ((NEWTON_GRAVITY_CGS*const_unit_mass_in_cgs*const_unit_time_in_cgs*const_unit_time_in_cgs/(const_unit_length_in_cgs*const_unit_length_in_cgs*const_unit_length_in_cgs)))
print 'G=', const_G
utherm                 = math.pi * const_G * surface_density * scale_height / (gamma-1)
v_disp                 = numpy.sqrt(2 * utherm)
soundspeed             = numpy.sqrt(utherm / (gamma * (gamma-1.)))
t_dyn                  = numpy.sqrt(scale_height / (const_G * surface_density))
t_cross                = scale_height / soundspeed
print 'dynamical time = ',t_dyn,' sound crossing time = ',t_cross,' sound speed= ',soundspeed,' 3D velocity dispersion = ',v_disp,' thermal_energy= ',utherm


# Parameters
periodic= 1            # 1 For periodic box
boxSize = 400.         #  [kpc]
boxSize_z = 800.
Radius  = 100.         # maximum radius of particles [kpc]
G       = const_G 

# File
fileName = "Disc-Patch.hdf5" 

#---------------------------------------------------
mass           = 1

#--------------------------------------------------


# using glass ICs
# read glass file and generate gas positions and tile it twice along z
inglass = 'glassCube_32.hdf5'
infile  = h5py.File(inglass, "r")
glass_p = infile["/PartType0/Coordinates"][:,:]
glass_h = infile["/PartType0/SmoothingLength"][:]

glass_p[:,:]      *= boxSize
glass_p[:,2]      -= boxSize_z/2
ndens_glass       = (glass_h.shape[0]) / (boxSize * boxSize * boxSize_z)
h_glass           = numpy.amin(glass_h) * (boxSize)

# random shuffle of glas ICs
numpy.random.seed(12345)
indx   = numpy.random.rand(numpy.shape(glass_h)[0])
indx   = numpy.argsort(indx)
glass_p = glass_p[indx, :]
glass_h = glass_h[indx]

# select numGas of them
numGas = 8192
pos    = glass_p[0:numGas,:]
pos2   = numpy.copy(pos)
pos2[:,2] += boxSize_z/2
pos = numpy.concatenate((pos, pos2), axis=0)
h      = glass_h[0:numGas]
h = numpy.append(h, h)
numGas = numpy.shape(pos)[0]

# compute further properties of ICs
column_density = fgas * surface_density * numpy.tanh(boxSize_z/2./scale_height)
enclosed_mass  = column_density * boxSize * boxSize
pmass          = enclosed_mass / numGas
meanrho        = enclosed_mass / boxSize*boxSize * boxSize_z
print 'pmass= ',pmass,' mean(rho) = ', meanrho,' entropy= ', (gamma-1) * utherm / meanrho**(gamma-1)

# desired density
rho            = surface_density / (2. * scale_height) / numpy.cosh(abs(pos[:,2])/scale_height)**2
u              = (1. + 0 * h) * utherm 
entropy        = (gamma-1) * u / rho**(gamma-1)
mass           = 0.*h + pmass
entropy_flag   = 0
vel            = 0 + 0 * pos

# Create and write output file

#File
file = h5py.File(fileName, 'w')

#Units
grp = file.create_group("/Units")
grp.attrs["Unit length in cgs (U_L)"] = const_unit_length_in_cgs
grp.attrs["Unit mass in cgs (U_M)"] = const_unit_mass_in_cgs 
grp.attrs["Unit time in cgs (U_t)"] = const_unit_length_in_cgs / const_unit_velocity_in_cgs
grp.attrs["Unit current in cgs (U_I)"] = 1.
grp.attrs["Unit temperature in cgs (U_T)"] = 1.

# Header
grp = file.create_group("/Header")
grp.attrs["BoxSize"] = [boxSize, boxSize, boxSize_z]
grp.attrs["NumPart_Total"] =  [numGas, 0, 0, 0, 0, 0]
grp.attrs["NumPart_Total_HighWord"] = [0, 0, 0, 0, 0, 0]
grp.attrs["NumPart_ThisFile"] = [numGas, 0, 0, 0, 0, 0]
grp.attrs["Time"] = 0.0
grp.attrs["NumFilesPerSnapshot"] = 1
grp.attrs["MassTable"] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
grp.attrs["Flag_Entropy_ICs"] = [entropy_flag]
grp.attrs["Dimension"] = 3

#Runtime parameters
grp = file.create_group("/RuntimePars")
grp.attrs["PeriodicBoundariesOn"] = periodic


# write gas particles
grp0   = file.create_group("/PartType0")

ds     = grp0.create_dataset('Coordinates', (numGas, 3), 'f')
ds[()] = pos

ds     = grp0.create_dataset('Velocities', (numGas, 3), 'f')
ds[()] = vel

ds     = grp0.create_dataset('Masses', (numGas,), 'f')
ds[()] = mass

ds     = grp0.create_dataset('SmoothingLength', (numGas,), 'f')
ds[()] = h

ds = grp0.create_dataset('InternalEnergy', (numGas,), 'f')
u = numpy.full((numGas, ), utherm)
if (entropy_flag == 1):
    ds[()] = entropy
else:
    ds[()] = u    

ids = 1 + numpy.linspace(0, numGas, numGas, endpoint=False)
ds = grp0.create_dataset('ParticleIDs', (numGas, ), 'L')
ds[()] = ids

print "Internal energy:", u[0]
file.close()
