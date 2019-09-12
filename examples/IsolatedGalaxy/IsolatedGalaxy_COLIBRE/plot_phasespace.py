import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import gridspec
import h5py
import sys
import glob
import numpy as np

import matplotlib as mpl

from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.patches as patches

tmin = 0.
tmax = 100.

red_patch   = patches.Patch(color='#ff0000', label='ionized')
green_patch = patches.Patch(color='#00ff00', label='neutral')
blue_patch  = patches.Patch(color='#0000ff', label='molecular')

logZmin =-2.
logZmax = 2.

norm = mpl.colors.Normalize(vmin=logZmin, vmax=logZmax)

def colorbar(mappable):
    ax = mappable.axes
    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    return fig.colorbar(mappable, cax=cax)

cmap = plt.cm.get_cmap('hsv')
cmap2 = plt.cm.get_cmap('coolwarm')

mH = 1.6726219e-24
Myr = 1.e6 * 365.25 * 24. * 3600.
alpha_B = 2.6e-13
Msol = 1.9891E33
kpc  = 3.086e21

## First, get list of all snapshots
reg_exp = "output_*.hdf5"
snaplist = sorted(glob.glob(reg_exp))

lognHmin = -6.
lognHmax =  8.
logTmin  =  1.
logTmax  =  9.5

nn = 100
vmin = -2.
vmax =  0. 

# for colorscale interpolation
varr = np.linspace(vmin, vmax, nn, endpoint = True)
yarr = np.linspace(0., 1., nn, endpoint = True)

isnap = 0
for snap in snaplist:
        fig = plt.figure()
        fig.subplots_adjust(bottom = 0.2, top = 0.8)

        fig.set_size_inches(6,3.)
        gs = gridspec.GridSpec(1,3, width_ratios = [0.9,0.9,0.1], wspace = 0.2)

        ax6 = plt.subplot(gs[0])
        ax6.set_title('Hydro variables')
        ax6.set_xlim(lognHmin, lognHmax)
        ax6.set_ylim(logTmin, logTmax)
        ax6.set_xlabel('log nH [cm$^{-3}$]')
        ax6.set_ylabel('log T [K]')

        ax8 = plt.subplot(gs[1])
        ax8.set_title('Subgrid variables')
        ax8.set_xlim(lognHmin, lognHmax)
        ax8.set_ylim(logTmin, logTmax)
        ax8.set_xlabel('log nH [cm$^{-3}$]')

        with h5py.File(snap, 'r') as f:
                dens        = f['PartType0/Densities'][:]
                dens_to_cgs = f['PartType0/Densities'].attrs['Conversion factor to physical CGS (including cosmological corrections)']
                temp        = f['PartType0/Temperatures'][:]
                massfracs   = f['PartType0/ElementMassFractions'][:]

                XH    = massfracs[:,0]
                lognH = np.log10(XH * dens * dens_to_cgs / mH)
                logT  = np.log10(temp)

                dens_subgrid = f['PartType0/SubgridDensity'][:]
                temp_subgrid = f['PartType0/SubgridTemperature'][:]

                lognH_sub = np.log10(XH * dens_subgrid * dens_to_cgs / mH)
                logT_sub  = np.log10(temp_subgrid)

                HIfrac  = f['PartType0/HydrogenNeutralFraction'][:]
                HIIfrac = f['PartType0/HydrogenIonizedFraction'][:]
                H2frac  = f['PartType0/HydrogenMolecularFraction'][:]

                logHI  = np.log10(HIfrac)
                logHII = np.log10(HIIfrac)
                logH2  = np.log10(2. * H2frac)

                time = f['Header'].attrs['Time']
                time_in_s = f['Units'].attrs['Unit time in cgs (U_t)']
                mass_in_g = f['Units'].attrs['Unit mass in cgs (U_M)']

                red   = np.interp(logHII, varr, yarr)
                green = np.interp(logHI , varr, yarr)
                blue  = np.interp(logH2 , varr, yarr)

                C = np.zeros((len(logHI),3))
                #C[:,0] = red[:]
                #C[:,1] = green[:]
                #C[:,2] = blue[:]
                C[:,0] = HIIfrac[:]
                C[:,1] = HIfrac[:]
                C[:,2] = 2. * H2frac[:]

        fig.suptitle('t = %6.1f Myr'%(time * time_in_s / Myr))

        ax6.scatter(lognH, logT, s = 2)
        ax8.scatter(lognH_sub, logT_sub, s = 2, c = C)

        ax8.legend(handles=[red_patch, green_patch, blue_patch], loc = 'upper right')
        outfile = 'plots/subgrid_%4.4i.png'%(isnap)
        print ('Plotting: %s'%(outfile))
        fig.savefig(outfile, dpi = 150)
        isnap = isnap + 1

        plt.close('all')
 





      
