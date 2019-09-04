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

red_patch = patches.Patch(color='red', label='HII stars')
grey_patch = patches.Patch(color='grey', label='HII gas')
black_patch = patches.Patch(color='black', label='expect.')

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

isnap = 0
for snap in snaplist:
        fig = plt.figure()
        fig.subplots_adjust(right = 0.85)

        fig.set_size_inches(6,2.5)
        gs = gridspec.GridSpec(1,3, width_ratios = [0.9,0.9,0.1], wspace = 0.1)

        ax6 = plt.subplot(gs[0])
        ax6.clear()
        ax6.set_xlim(-15,15)
        ax6.set_ylim(-15,15)
        ax6.set(adjustable='box-forced', aspect='equal')
        ax6.tick_params(labelleft = False, labelbottom = False)
        ax6.set_xticks([])
        ax6.set_yticks([])

        ax8 = plt.subplot(gs[1])
        ax8.set_xlim(-15,15)
        ax8.set_ylim(-15,15)
        ax8.set(adjustable='box-forced', aspect='equal')
        ax8.tick_params(labelleft = False, labelbottom = False)
        ax8.set_xticks([])
        ax8.set_yticks([])

        with h5py.File(snap, 'r') as f:
                StarID = f['PartType4/ParticleIDs'][:]
                HIImass = f['PartType4/HIIregions_mass_to_ionize'][:]
                tHII_Myr = f['PartType4/HIIregions_last_rebuild'][:]
                dens    = f['PartType4/BirthDensities'][:]
                mass    = f['PartType4/InitialMasses'][:]
                HIIstart = f['PartType4/BirthTimes'][:]
                dtstar   = f['PartType4/Timestep'][:]
                HIIgas   = f['PartType0/HIIregionsEndTime'][:]
                mass_gas = f['PartType0/Masses'][:]
                pos_gas_SU  = f['PartType0/Coordinates'][:]
                pos_stars_SU  = f['PartType4/Coordinates'][:]
                pos_conv = f['PartType0/Coordinates'].attrs['Conversion factor to physical CGS (including cosmological corrections)'] 

                dens_to_cgs = f['PartType4/BirthDensities'].attrs['Conversion factor to physical CGS (including cosmological corrections)'][0]
                XH = float( f['Parameters'].attrs['COLIBREChemistry:init_abundance_Hydrogen'] )
                Qbar = float( f['Parameters'].attrs['COLIBREFeedback:HIIregion_const_ionrate'] )
                dt_Myr = float (  f['Parameters'].attrs['COLIBREFeedback:HIIregion_rebuild_dt_Myr'] )
                HIIage_max_Myr = float ( f['Parameters'].attrs['COLIBREFeedback:HIIregion_maxage_Myr'] )

                time = f['Header'].attrs['Time']
                time_in_s = f['Units'].attrs['Unit time in cgs (U_t)']
                mass_in_g = f['Units'].attrs['Unit mass in cgs (U_M)']
                dtstar_in_Myr = dtstar * time_in_s / Myr

                StarAge_in_Myr = (time - HIIstart) * time_in_s / Myr

        fig.suptitle('t = %6.1f Myr'%(time * time_in_s / Myr))

        if isnap == 0:
                labstar = 'HII mass stars'
                labgas = 'HII mass gas'
                labexp = 'Expected'
        else:
                labstar = ''
                labgas = ''
                labexp = ''


        indxHIIstar   = np.where( (HIIstart >= 0) & (StarAge_in_Myr <= HIIage_max_Myr) )[0]
        HIImass_stars =  np.sum(HIImass[indxHIIstar])

        indxHIIgas  = np.where(HIIgas > 0.)[0]
        HIImass_gas = np.sum(mass_gas[indxHIIgas]) 

        cenx      = np.median(pos_gas_SU[:,0]) * pos_conv / kpc
        ceny      = np.median(pos_gas_SU[:,0]) * pos_conv / kpc
        pos_gas   = pos_gas_SU * pos_conv / kpc
        pos_stars = pos_stars_SU * pos_conv / kpc 


        ax6.scatter(pos_gas[indxHIIgas,0] - cenx, pos_gas[indxHIIgas,1]- ceny, s = 2, c = 'grey',  edgecolors = 'face')
        ax6.scatter(pos_stars[indxHIIstar,0] - cenx, pos_stars[indxHIIstar,1] - ceny, s = 4, c = 'red', marker = '*',  edgecolors = 'face')

        img8 = ax8.scatter(pos_stars[indxHIIstar,0] - cenx, pos_stars[indxHIIstar,1] - ceny, s = 8, marker = '*', cmap = cmap2, vmin = -2, vmax = 2, edgecolors = 'face', c = np.log10(HIImass[indxHIIstar] / mass[indxHIIstar]))

        #if (len (indxHIIstar) > 0):
        #        cbar = plt.colorbar(img8, ax = ax8)
        #        cbar.set_label('Nr neighbors')

        age_in_s = (time - HIIstart) * time_in_s
        t_half   = age_in_s + 0.5 * dt_Myr * Myr
        nH       = dens * dens_to_cgs / mH * XH
                
        indxnew = np.where(HIIstart >= 0)[0]

        HIImass_exp = np.zeros_like(HIImass[indxnew])

        for i in range (len(indxnew)):
                HIImass_exp[i] = 0.84 * mass[indxnew[i]] * (1. - np.exp( - alpha_B * nH[indxnew[i]] * (tHII_Myr[indxnew[i]] + 0.5 * dt_Myr) * Myr) ) * \
                ( 10. / nH[indxnew[i]] ) * (Qbar / 1.e12)
                if (HIImass[indxnew[i]] / HIImass_exp[i] < 0.9 or HIImass[indxnew[i]] / HIImass_exp[i] > 1.1):
                        print ('Problem: snapshot = %i\t HIImass/HIImass_exp = %.4f\t HIImass = %.4e\t HIImass_exp = %.4e\t Age [Myr] = %.4e\t%i\t%i'%(isnap,\
                                HIImass[indxnew[i]] / HIImass_exp[i], HIImass[indxnew[i]], HIImass_exp[i], \
                                age_in_s[indxnew[i]] / Myr, \
                                indxnew[i], StarID[indxnew[i]]))
                        print ('     More info: Init star mass = %.4e\t nH = %.4e\t HII last rebuild [Myr] = %.4f\t Star age = %.4f'%(mass[indxnew[i]], \
                               nH[indxnew[i]], tHII_Myr[indxnew[i]] + 0.5 * dt_Myr, age_in_s[indxnew[i]]/Myr))

        ax = plt.subplot(gs[2])
        cb1 = mpl.colorbar.ColorbarBase(ax, cmap=cmap2,norm=norm)
        cb1.set_label ('log nr neigh')
        
        outfile = 'plots/sml_HIIregion_%4.4i.png'%(isnap)
        print ('Plotting: %s'%(outfile))
        fig.savefig(outfile, dpi = 150)
        isnap = isnap + 1

        plt.close('all')
 
















      
