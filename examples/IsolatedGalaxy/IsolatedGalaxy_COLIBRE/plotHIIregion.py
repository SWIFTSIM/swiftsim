import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import gridspec
import h5py
import sys
import glob
import numpy as np

from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.patches as patches


def colorbar(mappable):
    ax = mappable.axes
    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    return fig.colorbar(mappable, cax=cax)

cmap = plt.cm.get_cmap('hsv')

mH = 1.6726219e-24
Myr = 1.e6 * 365.25 * 24. * 3600.
alpha_B = 2.6e-13
Msol = 1.9891E33

## First, get list of all snapshots
reg_exp = "output_*.hdf5"
snaplist = sorted(glob.glob(reg_exp))

fig = plt.figure()
fig.suptitle("HII region tests")
fig.set_size_inches(8,9.5)
gs = gridspec.GridSpec(3,2, wspace = 0.35)
ax = plt.subplot(gs[0])
ax.set_xlabel('Age of star particle [Myr]')
ax.set_ylabel('Mass of HII region [SU] over expected mass')
ax.set_xlim(-5,55)
ax.set_ylim(-0.1, 1.1)
ax.set_title ('Check star particle properties')

ax2 = plt.subplot(gs[1])
ax2.set_xlabel('Simulation time [Myr]')
ax2.set_ylabel('Ionized mass, mass to ionize')
ax2.set_yscale('log')
ax2.set_title ('Check gas particle properties')

ax4 = plt.subplot(gs[2])
ax4.set_xlabel('Age of star particle [Myr]')
ax4.set_ylabel('Log mass [Msol]')
ax4.set_yscale('log')
ax4.set_xlim(-0.2,5)

ax3 = plt.subplot(gs[3])
ax3.set_xlabel('Simulation time [Myr]')
ax3.set_ylabel('Ionized mass / mass to ionize')

ax5 = plt.subplot(gs[4])
ax5.set_xlabel('Age of star particle [Myr]')
ax5.set_ylabel('Timestep [Myr]')
ax5.set_yscale('log')


first = True
first2 = True
for snap in snaplist[0:-1]:
        with h5py.File(snap, 'r') as f:
                HIImass = f['PartType4/HIIregions_mass_to_ionize'][:]
                dens    = f['PartType4/BirthDensities'][:]
                mass    = f['PartType4/InitialMasses'][:]
                HIIstart = f['PartType4/BirthTimes'][:]
                dtstar   = f['PartType4/Timestep'][:]
                HIIgas   = f['PartType0/HIIregionsEndTime'][:]
                mass_gas = f['PartType0/Masses'][:]

                dens_to_cgs = f['PartType4/BirthDensities'].attrs['Conversion factor to physical CGS (including cosmological corrections)'][0]
                XH = float( f['Parameters'].attrs['COLIBREChemistry:init_abundance_Hydrogen'] )
                Qbar = float( f['Parameters'].attrs['COLIBREFeedback:HIIregion_const_ionrate'] )
                dt_Myr = float (  f['Parameters'].attrs['COLIBREFeedback:HIIregion_rebuild_dt_Myr'])

                time = f['Header'].attrs['Time']
                time_in_s = f['Units'].attrs['Unit time in cgs (U_t)']
                mass_in_g = f['Units'].attrs['Unit mass in cgs (U_M)']
                dtstar_in_Myr = dtstar * time_in_s / Myr


        age_in_s = (time - HIIstart) * time_in_s
        t_half   = age_in_s + 0.5 * dt_Myr * Myr
        nH       = dens * dens_to_cgs / mH * XH
                
        indxnew = np.where(HIIstart >= 0)[0]

        HIImass_exp = np.zeros_like(HIImass[indxnew])

        for i in range (len(indxnew)):
                HIImass_exp[i] = 0.84 * mass[indxnew[i]] * (1. - np.exp( - alpha_B * nH[indxnew[i]] * t_half[indxnew[i]] ) ) * \
                ( 10. / nH[indxnew[i]] ) * (Qbar / 1.e12)

        img1 = ax.scatter(age_in_s[indxnew] / Myr, HIImass[indxnew]/HIImass_exp, label = '%s'%(snap))
        ax.text(10, 0.7, 'shown in plot below', color = 'grey') 
        ax.add_patch( patches.Rectangle((-0.2, -0.1), 5.2, 1.05, fill=False, color = 'grey', linestyle = 'dotted' ) ) 

        indxprob = np.where( HIImass[indxnew]/HIImass_exp < 0.95 )[0] 
        if (len(indxprob) > 0):
                if first2:
                        ax4.scatter(age_in_s[indxnew[indxprob]] / Myr, HIImass[indxnew[indxprob]] * mass_in_g / Msol, label = 'Gas mass to ionize', color = 'black')
                        ax4.scatter(age_in_s[indxnew[indxprob]] / Myr, HIImass_exp[indxprob] * mass_in_g / Msol, label = 'Expected gas mass (ref)', color = 'red')
                        first2 = False
                else:
                        ax4.scatter(age_in_s[indxnew[indxprob]] / Myr, HIImass[indxnew[indxprob]] * mass_in_g / Msol, color = 'black')
                        ax4.scatter(age_in_s[indxnew[indxprob]] / Myr, HIImass_exp[indxprob] * mass_in_g / Msol, color = 'red')

        indxgas = np.where(HIIgas > 0.)[0]
        if (len(indxgas) > 0):
                ax3.axhline(1., color  = 'grey', linestyle = 'dotted')
                img2 = ax3.scatter(time * time_in_s / Myr, np.sum(mass_gas[indxgas]) / np.sum(mass[indxnew]), c = np.log10(len(indxgas)), vmin = 0, vmax = 4, cmap = cmap)
                if first:
                        ax2.scatter(time * time_in_s / Myr, np.sum(mass_gas[indxgas]), c = 'black', label = 'Ionized gas mass')
                        ax2.scatter(time * time_in_s / Myr, np.sum(mass[indxnew]), c = 'red', label = 'Gas mass to ionize (ref)')
                        first = False
                else:
                        ax2.scatter(time * time_in_s / Myr, np.sum(mass_gas[indxgas]), c = 'black')
                        ax2.scatter(time * time_in_s / Myr, np.sum(mass[indxnew]), c = 'red')


        ax5.scatter(age_in_s[indxnew] / Myr, dtstar_in_Myr[indxnew], color = 'black')

        if (len(indxnew) > 0):
                print('%s: min = %.4e, max = %.4e'%(snap, (HIImass[indxnew]/HIImass_exp).min(), (HIImass[indxnew]/HIImass_exp).max() ))


cbar = colorbar(img2)
cbar.set_label('log nr of particles')
ax2.legend()
ax4.legend()

fig.savefig('HIIregion_test.png', dpi = 150)

                                
