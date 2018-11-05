#!/usr/bin/env python
import numpy as np
import h5py 
import matplotlib.pyplot as plt
from scipy.integrate import odeint


lengthrun = 2001
numbpar = 5

radius = np.zeros((numbpar,lengthrun))
time = np.zeros(lengthrun)
for i in range(0,lengthrun):
    Data = h5py.File('hernquist_%04d.hdf5'%i,'r')
    header = Data['Header']
    time[i] = header.attrs['Time']
    particles = Data['PartType1']
    positions = particles['Coordinates']
    radius[:,i] = positions[:,0]-200.

col = ['b','r','c','y','k']

for i in range(0,numbpar):
    plt.plot(time,radius[i,:],col[i])
    plt.axhline(np.max(radius[i,:]),color=col[i],linestyle='--')
    plt.axhline(-np.max(radius[i,:]),color=col[i],linestyle='--')


plt.ylabel('Radial distance (kpc)')
plt.xlabel('Simulation time (internal units)')
plt.savefig('radial_infall.png')
plt.close()


