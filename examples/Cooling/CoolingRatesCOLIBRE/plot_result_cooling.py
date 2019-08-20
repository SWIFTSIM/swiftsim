
import sys
#---------------------------------------------------------------
if len(sys.argv) > 1:
	filename = sys.argv[1]
	z = float(sys.argv[2])
	d = float(sys.argv[3])
else:
	filename = 'cooling_output_lognH-1.00.dat'
	z = 0.
	d = -1.

outputfile = 'cooling_lognH%.2f_z%.2f.png'%(d, z)
#---------------------------------------------------------------
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import gridspec
import os

import h5py
import numpy as np
import warnings
warnings.filterwarnings('ignore')

SMALL_SIZE = 12
VERYSMALL_SIZE = SMALL_SIZE - 3
MEDIUM_SIZE = SMALL_SIZE + 2
BIG_SIZE    = SMALL_SIZE + 4
BIGGER_SIZE = SMALL_SIZE + 6
plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=VERYSMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIG_SIZE)  # fontsize of the figure title

cmap = plt.cm.get_cmap('viridis')
xmin = 1.
xmax = 9.5
ymin = -35.
ymax = -18.
xlab = 'log T [K]'

#newtables
primcool = [0,1,12,14,15,17,18]
atomcool = [2,3,4,5,6,7,8,9,10]
restcool = [11,13,16,19]

#################
# hardcoded !!!!
tables = '/cosma7/data/dp004/dc-ploe1/CoolingTables/2019_04/UV_dust1_CR1_G1_shield1.hdf5'
runname = 'UV_dust1_CR1_G1_shield1'
titlestring = '%s, z=%4.1f, log n [cm$^{-3}$]=%4.1f, log Z/Z$_{\odot}$=%4.1f'%(runname, float(z), float(d), 0.)
#################


with h5py.File(tables, "r") as f:
	IdentifierCooling = f['IdentifierCooling'].value

data = np.loadtxt(filename)

xx = np.log10(data[:,0])
Cool1D = np.log10(-1. * (data[:,2:]))
Cool1D_total = np.log10(-1. * (data[:,1]))

Cool1D_totalprim = np.zeros_like(Cool1D[:,0])
for i in range(len(primcool)):
	Cool1D_totalprim[:] = Cool1D_totalprim[:] + np.power(10., Cool1D[:,primcool[i]])
Cool1D_totalprim = np.log10(Cool1D_totalprim)

Cool1D_totalmetal = np.zeros_like(Cool1D[:,0])
for i in range(len(atomcool)):
	Cool1D_totalmetal[:] = Cool1D_totalmetal[:] + np.power(10., Cool1D[:,atomcool[i]])
for i in range(len(restcool)):
	Cool1D_totalmetal[:] = Cool1D_totalmetal[:] + np.power(10., Cool1D[:,restcool[i]])
Cool1D_totalmetal = np.log10(Cool1D_totalmetal)

fig = plt.figure()
fig.set_size_inches(10,6.2,forward=True)
fig.suptitle(titlestring)
fig.subplots_adjust(left = 0.1, right = 0.95, bottom = 0.15, top = 0.7)
gs = gridspec.GridSpec(1,3,wspace=0, hspace=0)

# Cooling plot
ax = plt.subplot(gs[0])
ax.set_xlim(xmin,xmax)
ax.set_ylim(ymin,ymax)
ax.set_ylabel('log $\Lambda$/$n_{\mathrm{H}}^2$ [erg cm$^{3}$ s$^{-1}$]')

ax.set_xlabel(xlab)
ax.xaxis.set_ticks(np.arange(xmin, xmax, 1.))
ax.plot(xx, Cool1D_total, color = 'grey', lw = 2, label = 'Total')

icount = 0
c = cmap(np.linspace(0,1,len(primcool)))
for i in range (len(primcool)):
    icool = primcool[i]
    linec = c[icount]
    if icount%3 == 0:
       lines = '-'
    elif icount%3 == 1:
        lines = '--'
    else:
        lines = '-.'

    ax.plot(xx, Cool1D[..., icool], color = linec, lw = 2, ls = lines, label = '%s'%(IdentifierCooling[icool].decode('utf-8')))
    icount = icount + 1

ax.plot(xx, Cool1D_totalprim, color = 'black', lw = 2, ls = ':', label = 'TotalPrim')

handles, labels = ax.get_legend_handles_labels()
ax.legend(handles, labels, bbox_to_anchor=(0., 1.02, 1., .102), loc=3, mode = 'expand', borderaxespad=0., ncol = 2, fontsize = VERYSMALL_SIZE, handlelength = 4)


ax = plt.subplot(gs[1])
ax.set_xlim(xmin,xmax)
ax.set_ylim(ymin,ymax)
ax.tick_params(labelleft='off')

ax.set_xlabel(xlab)
ax.xaxis.set_ticks(np.arange(xmin, xmax, 1.))
ax.plot(xx, Cool1D_total, color = 'grey', lw = 2, label = 'Total')
icount = 0
c = cmap(np.linspace(0,1,len(atomcool)))
for i in range (len(atomcool)):
    icool = atomcool[i]
    linec = c[icount]
    if icount%3 == 0:
        lines = '-'
    elif icount%3 == 1:
        lines = '--'
    else:
        lines = '-.'

    ax.plot(xx, Cool1D[..., icool], color = linec, lw = 2, ls = lines, label = '%s'%(IdentifierCooling[icool].decode('utf-8')))
    icount = icount + 1

ax.plot(xx, Cool1D_totalmetal, color = 'black', lw = 2, ls = ':', label = 'TotalMetal')

handles, labels = ax.get_legend_handles_labels()
ax.legend(handles, labels, bbox_to_anchor=(0., 1.02, 1., .102), loc=3, mode = 'expand', borderaxespad=0., ncol = 2, fontsize = VERYSMALL_SIZE, handlelength = 4)


ax = plt.subplot(gs[2])
ax.set_xlim(xmin,xmax)
ax.set_ylim(ymin,ymax)
ax.tick_params(labelleft='off')

ax.set_xlabel(xlab)
ax.xaxis.set_ticks(np.arange(xmin, xmax, 1.))
ax.plot(xx, Cool1D_total, color = 'grey', lw = 2, label = 'Total')
icount = 0
c = cmap(np.linspace(0,1,len(restcool)))
for i in range (len(restcool)):
    icool = restcool[i]
    linec = c[icount]
    if icount%3 == 0:
        lines = '-'
    elif icount%3 == 1:
        lines = '--'
    else:
        lines = '-.'

    ax.plot(xx, Cool1D[..., icool], color = linec, lw = 2, ls = lines, label = '%s'%(IdentifierCooling[icool].decode('utf-8')))
    icount = icount + 1

ax.plot(xx, Cool1D_totalmetal, color = 'black', lw = 2, ls = ':', label = 'TotalMetal')

handles, labels = ax.get_legend_handles_labels()
ax.legend(handles, labels, bbox_to_anchor=(0., 1.02, 1., .102), loc=3, mode = 'expand', borderaxespad=0., ncol = 2, fontsize = VERYSMALL_SIZE, handlelength = 4)

fig.savefig(outputfile, dpi = 100)
print ('Figure saved as: %s'%(outputfile))
plt.close('all')

