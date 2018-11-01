from galpy.potential import NFWPotential
from galpy.orbit import Orbit
from galpy.util import bovy_conversion
import numpy as np
import matplotlib.pyplot as plt
from astropy import units
import write_gadget as wg
import h5py as h5

C = 8.0
M_200 = 2.0
N_PARTICLES = 1
def write_hdf5():
	print("\nInitial conditions written to 'test_nfw.hdf5'")

	pos=np.array([8.0,0.0,0.0])+500.0
	vel=np.array([0.0, 240.0, 5.0])
	ids=np.array([1.0])
	mass=np.array([1.0])

    #File
	file = h5.File("test_nfw.hdf5", 'w')

	#Units
	grp = file.create_group("/Units")
	grp.attrs["Unit length in cgs (U_L)"] = 3.086e+21
	grp.attrs["Unit mass in cgs (U_M)"] = 1.988e+33 
	grp.attrs["Unit time in cgs (U_t)"] = 3.086e+16
	grp.attrs["Unit current in cgs (U_I)"] = 1.
	grp.attrs["Unit temperature in cgs (U_T)"] = 1.

	# Header
	grp = file.create_group("/Header")
	grp.attrs["BoxSize"] = 1000.0
	grp.attrs["NumPart_Total"] =  [0, N_PARTICLES, 0, 0, 0, 0]
	grp.attrs["NumPart_Total_HighWord"] = [0, 0, 0, 0, 0, 0]
	grp.attrs["NumPart_ThisFile"] = [0, N_PARTICLES, 0, 0, 0, 0]
	grp.attrs["Time"] = 0.0
	grp.attrs["NumFilesPerSnapshot"] = 1
	grp.attrs["MassTable"] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
	grp.attrs["Flag_Entropy_ICs"] = [0, 0, 0, 0, 0, 0]
	grp.attrs["Dimension"] = 3

	#Runtime parameters
	grp = file.create_group("/RuntimePars")
	grp.attrs["PeriodicBoundariesOn"] = 1

	#Particle group
	grp1 = file.create_group("/PartType1")
	ds = grp1.create_dataset('Velocities', (N_PARTICLES, 3), 'f')
	ds[()] = vel

	ds = grp1.create_dataset('Masses', (N_PARTICLES,), 'f')
	ds[()] = mass

	ds = grp1.create_dataset('ParticleIDs', (N_PARTICLES, ), 'L')
	ds[()] = ids

	ds = grp1.create_dataset('Coordinates', (N_PARTICLES, 3), 'd')
	ds[()] = pos

	file.close()	

write_hdf5()
