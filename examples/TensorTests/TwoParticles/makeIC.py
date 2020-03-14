import numpy as np
import h5py

npart = 2
pos = np.array([[1.,1.,1.], [2.,1.,1.]])
#pos = np.array([[1.,1.,1.], [1.57735,1.57735,1.57735]]) # diagonals of T are 0
vel = np.array([[0.,0.,0.], [0.,0.,0.]])
ids = np.arange(1,npart+1)
mass = np.array([1.]*npart)

# File
file = h5py.File("test_tensors.hdf5", "w")

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
