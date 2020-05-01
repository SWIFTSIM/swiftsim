# tidal_strengths.py
# plot the tidal field strengths from the tidal tensors

import os, sys, numpy as np, matplotlib as mpl, h5py
import matplotlib.pyplot as plt
from numpy import linalg as LA

SOLAR_MASS = 1.989e33
CM_PER_KPC = 3.085678e21
SEC_PER_GIGAYEAR = 3.155e16

dpi = 200

print '*** tidal strengths ***'
print

f = h5py.File( sys.argv[1], 'r' )
G = float(f['/Parameters'].attrs['PhysicalConstants:G'])
eps = float(f['/Parameters'].attrs['Gravity:max_physical_baryon_softening'])
mesh = float(f['/Parameters'].attrs['Gravity:mesh_side_length'])
BoxSize = np.array(f['/Header'].attrs['BoxSize'])[0]
print 'BoxSize', BoxSize
print 'mesh_side_length',mesh

U_M = float(f['/Units'].attrs['Unit mass in cgs (U_M)'])
U_L = float(f['/Units'].attrs['Unit length in cgs (U_L)'])
U_t = float(f['/Units'].attrs['Unit time in cgs (U_t)'])
print 'Unit mass',U_M/SOLAR_MASS,'Msun'
print 'Unit length',U_L/CM_PER_KPC,'kpc'
print 'Unit time',U_t/SEC_PER_GIGAYEAR,'Gyr'

pos = f['/PartType4/Coordinates'][:]
mass = f['/PartType4/Masses'][:]
tidal_tensor = f['/PartType4/TidalTensors'][:]
tt_to_cgs = f["PartType4/TidalTensors"].attrs[
    "Conversion factor to physical CGS (including cosmological corrections)"
][0]
tidal_tensor *= tt_to_cgs * (SEC_PER_GIGAYEAR)**2

f.close()

# Just get the last one
tidal_tensor = tidal_tensor[:,-6:]


Npart = len(pos)
print 'Number of particles:',Npart

pos[:,0] -= BoxSize/2.
pos[:,1] -= BoxSize/2.
pos[:,2] -= BoxSize/2.

# centre on the galaxy
print np.mean(pos[:,0]), np.mean(pos[:,1]), np.mean(pos[:,2])
pos[:,0] -= np.mean(pos[:,0])
pos[:,1] -= np.mean(pos[:,1])
pos[:,2] -= np.mean(pos[:,2])

r = np.linalg.norm(pos,axis=1)
print 'r min/max:', np.min(r), np.max(r)

# Get the tidal strengths (eigenvalues of tensor) for each particle
tidal_strength1 = []
tidal_strength2 = []
tidal_strength3 = []
for i in xrange(len(tidal_tensor)):
  T = np.zeros((3,3))
  T[0,0] = tidal_tensor[i][0]
  T[0,1] = tidal_tensor[i][1]
  T[0,2] = tidal_tensor[i][2]
  T[1,1] = tidal_tensor[i][3]
  T[1,2] = tidal_tensor[i][4]
  T[2,2] = tidal_tensor[i][5]
  T[1,0] = tidal_tensor[i][1]
  T[2,0] = tidal_tensor[i][2]
  T[2,1] = tidal_tensor[i][4]
  eigenval = LA.eigvalsh(T)
  tidal_strength1.append( max(eigenval) )
  tidal_strength3.append( min(eigenval) )
  if max(eigenval) == min(eigenval):
    tidal_strength2.append( 0. )
  else:
    tidal_strength2.append( max(eigenval[eigenval<max(eigenval)]) )
  #if r[i]>20: 
  #  w, v = LA.eigh(T)
  #  print i, r[i], v[:,0], v[:,1], v[:,2]
tidal_strength1 = np.array(tidal_strength1) 
tidal_strength2 = np.array(tidal_strength2) 
tidal_strength3 = np.array(tidal_strength3) 


# The Plummer sphere
M = 1.0e10 # Msun
Rc = 2. # kpc
rpl = np.linspace(0.0,30,500)
rho = 3.*M/(4.*np.pi*Rc**3) * (1.+rpl**2/Rc**2)**-2.5
a = -G*M*(rpl/Rc)/Rc**2*(1+(rpl/Rc)**2)**-1.5
Tpl_r = G*M/Rc**3 * (2.*(rpl/Rc)**2-1.) / (1+(rpl/Rc)**2)**2.5
Tpl_phi = -G*M/Rc**3 / (1+(rpl/Rc)**2)**1.5

# Renaud+11
lambda1_plus_Omega = G*M/Rc**3 / (1+(rpl/Rc)**2)**2.5 * 3.*(rpl/Rc)**2 

SHOW_OMEGA=True

plt.figure(figsize=(6,4.5))
if SHOW_OMEGA:
  plt.plot(r, tidal_strength1-0.5*(tidal_strength2+tidal_strength3), '.', c='cornflowerblue', ms=2, mec='none')
plt.plot(r, tidal_strength3, '.', c='gold', ms=2, mec='none')
plt.plot(r, tidal_strength2, '.', c='gold', ms=2, mec='none')
plt.plot(r, tidal_strength1, '.', c='orangered', ms=2, mec='none')
plt.plot(rpl,lambda1_plus_Omega, 'k-' , lw=2, label='$\partial^2\Phi_\mathrm{pl} / \partial r^2 + \Omega^2$')
plt.plot(rpl,  Tpl_r, 'k--', lw=2, label='$\partial^2\Phi_\mathrm{pl} / \partial r^2$')
plt.plot(rpl,Tpl_phi, 'k-.', lw=2, label='$\partial^2\Phi_\mathrm{pl} / \partial [\phi,\\theta]^2$')
plt.plot([], [], '.', c='orangered', ms=5, mec='none', label='$\lambda_{1}$')
plt.plot([], [], '.', c='gold', ms=5, mec='none', label='$\lambda_2$, $\lambda_3$')
if SHOW_OMEGA:
  plt.plot([], [], '.', c='cornflowerblue', ms=5, mec='none', label='$\lambda_1-0.5(\lambda_2+\lambda_3)$')
plt.legend(loc=4, numpoints=1)
plt.xlabel('r [kpc]')
plt.ylabel('Tidal field strength [Gyr$^{-2}$]')
plt.xscale('log')
plt.xlim(eps,30)
plt.ylim(-8000,4000)
plt.yscale('symlog', linthreshy=100)
plt.tight_layout(pad=0.5)

exp = np.floor(np.log10(Npart))
figname = 'plummer_N%dE%d_mesh%d_L%d.png' % (Npart/10**exp, exp, mesh, BoxSize)
print 'figname:',figname
plt.savefig(figname, dpi=dpi)

#plt.show()


