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

sdot = 4
Myr = 1.e6 * 365.25 * 24. * 3600.

## First, get list of all snapshots
reg_exp = "output_*.hdf5"
snaplist = sorted(glob.glob(reg_exp))

snapshottimes_Myr = []
ratioHIIgas_to_stars = []

for isnap in range( len(snaplist) ):
        snap = snaplist[isnap]


        with h5py.File(snap, 'r') as f:
                StarID = f['PartType4/ParticleIDs'][:]
                HIImass = f['PartType4/HIIregions_mass_to_ionize'][:]
                kernelmass = f['PartType4/HIIregions_mass_in_kernel'][:]
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

                Gas_StarID = f['PartType0/HIIregionsStarID'][:]

                dens_to_cgs = f['PartType4/BirthDensities'].attrs['Conversion factor to physical CGS (including cosmological corrections)'][0]
                XH = float( f['Parameters'].attrs['COLIBREChemistry:init_abundance_Hydrogen'] )
                dt_Myr = float (  f['Parameters'].attrs['COLIBREFeedback:HIIregion_rebuild_dt_Myr'] )
                HIIage_max_Myr = float ( f['Parameters'].attrs['COLIBREFeedback:HIIregion_maxage_Myr'] )

                time = f['Header'].attrs['Time']
                time_in_s = f['Units'].attrs['Unit time in cgs (U_t)']
                mass_in_g = f['Units'].attrs['Unit mass in cgs (U_M)']
                dtstar_in_Myr = dtstar * time_in_s / Myr

                StarAge_in_Myr = (time - HIIstart) * time_in_s / Myr


        snapshottimes_Myr.append(time * time_in_s / Myr)

        MedGasPartMass = np.median(mass_gas)

        indx_star = np.where(HIImass > 0.)[0]

        NrList      = []
        StarHIImass = []
        GasHIImass  = []
        ngbmass     = []

        IndxGasInHII = []


        # loop over all star particles that are producing an HII region right now
        for i in range( len(indx_star) ):
                NrList.append(i)
                # find gas particles that are in the HII region of this star
                indx_gas = np.where (Gas_StarID == StarID[indx_star[i]])[0]
                StarHIImass.append(HIImass[indx_star[i]])
                GasHIImass.append(np.sum(mass_gas[indx_gas]))
                ngbmass.append(kernelmass[indx_star[i]])
                #print (NrList[-1], StarHIImass[-1], GasHIImass[-1], ngbmass[-1])
                for ii in range (len(indx_gas)):
                      IndxGasInHII.append(indx_gas[ii])

        ratioHIIgas_to_stars.append(np.sum(GasHIImass)/np.sum(StarHIImass))
        
        ii = np.where(np.asarray(GasHIImass) <= 0.)[0]

        fig = plt.figure()
        fig.set_size_inches(8,7)
        fig.suptitle('t = %6.1f Myr, HII gas/stars = %.4f'%(time * time_in_s / Myr, np.sum(GasHIImass)/np.sum(StarHIImass)))
        gs = gridspec.GridSpec(2,2, wspace = 0.3, hspace = 0.35)

        ax = plt.subplot(gs[0])
        ax.set_title('Gas in HII regions over expectation')
        ax.set_xlabel('HII region number')
        ax.set_ylabel('log HII mass (gas) / HII mass (*)')
        ax.set_ylim(-2., 3.5)
        ax.scatter(NrList, np.log10( np.asarray(GasHIImass) / np.asarray(StarHIImass)), s = sdot)
        yy = np.zeros(len(ii))
        yy[:] = -1.5
        ax.scatter(np.asarray(NrList)[ii], yy, s = sdot)
        ax.axhline(0., linestyle = 'dotted')

        #ax = plt.subplot(gs[1])
        #ax.set_xlabel('HII region number')
        #ax.set_ylabel('HII mass (*)')
        #ax.scatter(np.asarray(NrList), np.asarray(StarHIImass))

        #ax = plt.subplot(gs[2])
        #ax.set_xlabel('HII region number')
        #ax.set_ylabel('Kernel mass (*)')
        #ax.scatter(NrList, ngbmass)

        ax = plt.subplot(gs[1])
        ax.set_title('Nr particles in kernel')
        ax.set_xlabel('HII region number')
        ax.set_ylabel('Kernel mass (*) / part mass')
        ax.set_yscale('log')
        ax.set_ylim(0.5, 300.)
        ax.scatter(NrList, ngbmass / MedGasPartMass, s = sdot) 
        ax.axhline(48, linestyle = 'dotted')

        ax = plt.subplot(gs[2])
        ax.set_title('Number of gas particles to heat')
        ax.set_xlabel('HII region number')
        ax.set_ylabel('HII mass (*) / part mass')
        ax.set_yscale('log')
        ax.axhline(1., linestyle = 'dotted')
        ax.set_ylim(1.e-4, 1000.)
        ax.scatter(NrList, np.asarray(StarHIImass) / MedGasPartMass, s = sdot)

        ax = plt.subplot(gs[3])
        ax.set_title('Number of gas particles heated')
        ax.set_xlabel('HII region number')
        ax.set_ylabel('HII mass (gas) / part mass')
        ax.set_ylim(1.e-4, 1000.)
        ax.set_yscale('log')
        ax.scatter(NrList, np.asarray(GasHIImass)  / MedGasPartMass, s = sdot)
        yy = np.zeros(len(ii))
        yy[:] = 5.e-4
        ax.scatter(np.asarray(NrList)[ii], yy, s = sdot)

        plt.subplots_adjust(top=0.9)
        fig.savefig('plots/IndivHIIregion_%4.4i.png'%(isnap), dpi = 150)
        print ('plots/IndivHIIregion_%4.4i.png'%(isnap))

        plt.close('all')



fig = plt.figure()
gs = gridspec.GridSpec(1,1)

ax = plt.subplot(gs[0])
ax.set_title('Gas mass heated / Gas mass to heat')
ax.set_xlabel('Simulation time [Myr]')
ax.set_ylabel('Log heated mass / to heat mass')
ax.scatter(snapshottimes_Myr, np.log10(np.asarray(ratioHIIgas_to_stars)), s = sdot)
ax.axhline(0.)


fig.savefig('plots/overview.png', dpi = 150)
print ('plots/overview.png')

      
