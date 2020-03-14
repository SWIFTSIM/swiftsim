import os
import numpy as np
import h5py

BoxSize = 1000.

print 'Making ICs'

data = np.loadtxt('plummer.dat', skiprows=1)
npart = len(data)
print 'Plummer sphere with %d particles'%(npart)

mass = data[:,1]           # 1e10 Msun
pos = data[:,2:5]/ 1e3     # kpc
vel = data[:,5:8]* 1.02273 # km/s
ids = np.arange(1,npart+1)

pos[:,0] += BoxSize/2.
pos[:,1] += BoxSize/2.
pos[:,2] += BoxSize/2.

# File
file = h5py.File("plummer.hdf5", "w")

# Units
grp = file.create_group("/Units")
grp.attrs["Unit length in cgs (U_L)"] = 3.086e21
grp.attrs["Unit mass in cgs (U_M)"] = 1.988e33
grp.attrs["Unit time in cgs (U_t)"] = 3.086e16
grp.attrs["Unit current in cgs (U_I)"] = 1.0
grp.attrs["Unit temperature in cgs (U_T)"] = 1.0

# Header
grp = file.create_group("/Header")
grp.attrs["BoxSize"] = 1000.0
grp.attrs["NumPart_Total"] = [0, 0, 0, 0, npart, 0]
grp.attrs["NumPart_Total_HighWord"] = [0, 0, 0, 0, 0, 0]
grp.attrs["NumPart_ThisFile"] = [0, 0, 0, 0, npart, 0]
grp.attrs["Time"] = 0.0
grp.attrs["NumFilesPerSnapshot"] = 1
grp.attrs["MassTable"] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
grp.attrs["Flag_Entropy_ICs"] = [0, 0, 0, 0, 0, 0]
grp.attrs["Dimension"] = 3

# Particle group
grp1 = file.create_group("/PartType4")
ds = grp1.create_dataset("Velocities", (npart, 3), "f", data=vel)

ds = grp1.create_dataset("Masses", (npart,), "f", data=mass)

ds = grp1.create_dataset("ParticleIDs", (npart,), "L", data=ids)

ds = grp1.create_dataset("Coordinates", (npart, 3), "d", data=pos)

file.close()
