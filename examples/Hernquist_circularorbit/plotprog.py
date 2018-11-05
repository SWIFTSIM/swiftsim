#!/usr/bin/env python
import numpy as np
import h5py 
import matplotlib.pyplot as plt
from scipy.integrate import odeint

t = np.linspace(0,40,1e5)
y0 = [0,10]
a = 30.0
G = 4.300927e-06
M = 1e15
GM = G*M


lengthrun = 2001
numbpar = 3

radius = np.zeros((numbpar,lengthrun))
xx = np.zeros((numbpar,lengthrun))
yy = np.zeros((numbpar,lengthrun))
zz = np.zeros((numbpar,lengthrun))
time = np.zeros(lengthrun)
for i in range(0,lengthrun):
    print(i)
    Data = h5py.File('output_%04d.hdf5'%i,'r')
    header = Data['Header']
    time[i] = header.attrs['Time']
    particles = Data['PartType1']
    positions = particles['Coordinates']
    xx[:,i] = positions[:,0]-500.
    yy[:,i] = positions[:,1]-500.
    zz[:,i] = positions[:,2]-500.

col = ['b','r','c','y','k']
print(np.shape(xx), np.shape(yy), np.shape(zz))

for i in range(0,numbpar):
    plt.plot(xx[i,:],yy[i,:],col[i])

plt.ylabel('y (kpc)')
plt.xlabel('x (kpc)')
plt.savefig('xyplot.png')
plt.close()


for i in range(0,numbpar):
    plt.plot(xx[i,:],zz[i,:],col[i])

plt.ylabel('z (kpc)')
plt.xlabel('x (kpc)')
plt.savefig('xzplot.png')
plt.close()

for i in range(0,numbpar):
    plt.plot(yy[i,:],zz[i,:],col[i])

plt.ylabel('z (kpc)')
plt.xlabel('y (kpc)')
plt.savefig('yzplot.png')
plt.close()
