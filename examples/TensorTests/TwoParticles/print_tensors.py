### Calculate tidal tensors with numerical differentiation for verification ###

import os, sys
import numpy as np
from numpy import linalg as LA
import h5py

GRAVITY = 4.29944e-06

def eval_gravity(i, mass, pos_i, pos, eps):
  eps2 = eps**2
  acc = np.zeros(3)
  for j in xrange(len(pos)):
    if j==i: continue
    dx = pos[j][0]-pos_i[0]
    dy = pos[j][1]-pos_i[1]
    dz = pos[j][2]-pos_i[2]
    r2 = dx**2 + dy**2 + dz**2
    if r2 >= 3**2*eps2:
      fac = GRAVITY*mass[j] / r2**1.5
    else:
      fac = GRAVITY*mass[j] / (r2+eps2)**1.5
    acc[0] += dx * fac
    acc[1] += dy * fac
    acc[2] += dz * fac
  return acc


if len(sys.argv) != 2:
  print 'usage: %s snapshot.hdf5' % (sys.argv[0])
  exit()

f = h5py.File( sys.argv[1], 'r' )
eps = float(f['/Parameters'].attrs['Gravity:max_physical_baryon_softening'])
ids = f['/PartType4/ParticleIDs'][:]
pos = f['/PartType4/Coordinates'][:] 
mass = f['/PartType4/Masses'][:] 
T = f['/PartType4/TidalTensors'][:] 
f.close()


N = len(pos)

print 'ids =', ids
print 'pos =', pos
print 'mass =', mass
print 'eps =', eps

dx = 1e-8 * eps

for i in xrange(N):
  print 'Particle ',i

  print '  TT in snapshot (upper):'
  print ' ', T[i][-6:] # Only print last one

  acc = eval_gravity(i, mass, pos[i], pos, eps)
  extraacc = np.zeros((3,3))
  extraacc[0] = eval_gravity(i, mass, pos[i]+np.array([1.,0.,0.])*dx, pos, eps)
  extraacc[1] = eval_gravity(i, mass, pos[i]+np.array([0.,1.,0.])*dx, pos, eps)
  extraacc[2] = eval_gravity(i, mass, pos[i]+np.array([0.,0.,1.])*dx, pos, eps)

  tide = np.zeros((3,3))
  tide[0] = (extraacc[:,0] - acc[0]) / (dx)
  tide[1] = (extraacc[:,1] - acc[1]) / (dx)
  tide[2] = (extraacc[:,2] - acc[2]) / (dx)
  print '  Expected tidal tensor:'
  print ' ', tide

  eigenval = LA.eigvalsh(tide)
  print '  eigenvalues:', eigenval
