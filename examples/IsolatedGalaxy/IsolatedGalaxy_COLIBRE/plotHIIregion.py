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

tmin = 0.
tmax = 100.

red_patch = patches.Patch(color='red', label='HII stars')
grey_patch = patches.Patch(color='grey', label='HII gas')
black_patch = patches.Patch(color='black', label='expect.')


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
for snap in snaplist[:-1]:
        fig = plt.figure()
        fig.set_size_inches(12,9)
        gs = gridspec.GridSpec(3,3, wspace = 0.45, hspace = 0.35)
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
        ax2.set_xlim(tmin,tmax)

        ax4 = plt.subplot(gs[2])
        ax4.set_xlabel('Age of star particle [Myr]')
        ax4.set_ylabel('Log mass [Msol]')
        ax4.set_yscale('log')
        ax4.set_xlim(-0.2,5)

        ax3 = plt.subplot(gs[3])
        ax3.set_xlabel('Simulation time [Myr]')
        ax3.set_ylabel('Ionized mass / mass to ionize')
        ax3.set_xlim(tmin,tmax)

        ax5 = plt.subplot(gs[4])
        ax5.set_xlabel('Age of star particle [Myr]')
        ax5.set_ylabel('Timestep [Myr]')
        ax5.set_yscale('log')
        ax5.set_xlim(tmin,tmax)

        ax6 = plt.subplot(gs[5])
        ax6.clear()
        ax6.set_xlabel('xpos')
        ax6.set_ylabel('ypos')
        ax6.set_xlim(-15,15)
        ax6.set_ylim(-15,15)

        ax7 = plt.subplot(gs[6])
        ax7.set_xlabel('Age of star particle [Myr]')
        ax7.set_ylabel('Age last rebuild HII region [Myr]')
        ax7.set_xlim(-5,35)
        ax7.set_ylim(-5,35)

        ax8 = plt.subplot(gs[8])
        ax8.set_xlabel('xpos')
        ax8.set_ylabel('ypos')
        ax8.set_xlim(-15,15)
        ax8.set_ylim(-15,15)

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

        if (len (indxHIIstar) > 0):
                cbar = plt.colorbar(img8, ax = ax8)
                cbar.set_label('Nr neighbors')

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
        

        img1 = ax.scatter(age_in_s[indxnew] / Myr, HIImass[indxnew]/HIImass_exp, label = '%s'%(snap), color = 'black', s = 2)
        ax.text(10, 0.7, 'shown in plot below', color = 'grey') 
        ax.add_patch( patches.Rectangle((-0.2, -0.1), 5.2, 1.05, fill=False, color = 'grey', linestyle = 'dotted' ) ) 

        indxprob = np.where( HIImass[indxnew]/HIImass_exp < 0.95 )[0] 
        if (len(indxprob) > 0):
                ax4.scatter(age_in_s[indxnew[indxprob]] / Myr, HIImass[indxnew[indxprob]] * mass_in_g / Msol, label = labstar, color = 'black', s = 4)
                ax4.scatter(age_in_s[indxnew[indxprob]] / Myr, HIImass_exp[indxprob] * mass_in_g / Msol, label = labexp, color = 'red', s = 4)

        if (len(indxHIIgas) > 0):
                ax3.axhline(1., color  = 'grey', linestyle = 'dotted')
                img2 = ax3.scatter(time * time_in_s / Myr, HIImass_gas/HIImass_stars, c = np.log10(len(indxHIIgas)), vmin = 0, vmax = 4, cmap = cmap, s = 8,  edgecolors = 'face')
                ax2.scatter(time * time_in_s / Myr, HIImass_gas , c = 'grey', label = labgas, s = 2)
                ax2.scatter(time * time_in_s / Myr, HIImass_stars, c = 'red', label = labstar, s = 2)
                cbar = plt.colorbar(img2, ax = ax3)
                cbar.set_label('log nr of particles')


        ax5.scatter(age_in_s[indxnew] / Myr, dtstar_in_Myr[indxnew], color = 'red', s= 2)

        if (len(indxnew) > 0):
                print('%s: min = %.4e, max = %.4e'%(snap, (HIImass[indxnew]/HIImass_exp).min(), (HIImass[indxnew]/HIImass_exp).max() ))

        ax.legend(handles=[red_patch, grey_patch, black_patch], loc = 'lower right')

        ax7.plot([-5,35], [-5,35], linestyle = 'dotted', color = 'black')
        ax7.scatter(age_in_s[indxnew] / Myr, tHII_Myr[indxnew], color = 'red', s = 4)

        fig.savefig('plots/HIIregion_%4.4i.png'%(isnap), dpi = 150)
        isnap = isnap + 1

        plt.close('all')
 
















      
