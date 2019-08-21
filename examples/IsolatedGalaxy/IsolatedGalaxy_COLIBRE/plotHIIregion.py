import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import h5py
import sys
import glob
import numpy as np

mH = 1.6726219e-24
Myr = 1.e6 * 365.25 * 24. * 3600.
alpha_B = 2.6e-13

## First, get list of all snapshots
reg_exp = "output_*.hdf5"
snaplist = sorted(glob.glob(reg_exp))

fig = plt.figure()
ax = plt.subplot(111)
ax.set_xlabel('Age of star particle [Myr]')
ax.set_ylabel('Mass of HII region [SU] over expected mass')

for snap in snaplist[0:-1]:
        with h5py.File(snap, 'r') as f:
                HIImass = f['PartType4/HIIregions_mass_to_ionize'][:]
                dens    = f['PartType4/BirthDensities'][:]
                mass    = f['PartType4/InitialMasses'][:]
                HIIstart = f['PartType4/BirthTimes'][:]

                dens_to_cgs = f['PartType4/BirthDensities'].attrs['Conversion factor to physical CGS (including cosmological corrections)'][0]
                XH = float( f['Parameters'].attrs['COLIBREChemistry:init_abundance_Hydrogen'] )
                Qbar = float( f['Parameters'].attrs['COLIBREFeedback:HIIregion_const_ionrate'] )

                time = f['Header'].attrs['Time']
                time_in_s = f['Units'].attrs['Unit time in cgs (U_t)']

        age_in_s = (time - HIIstart) * time_in_s
        nH       = dens * dens_to_cgs / mH * XH
                
        indxnew = np.where(HIIstart >= 0)[0]

        HIImass_exp = np.zeros_like(HIImass[indxnew])

        for i in range (len(indxnew)):
                HIImass_exp[i] = 0.84 * mass[indxnew[i]] * (1. - np.exp( - alpha_B * nH[indxnew[i]] * age_in_s[indxnew[i]] ) ) * \
                ( 10. / nH[indxnew[i]] ) * (Qbar / 1.e12)

        ax.scatter(age_in_s[indxnew] / Myr, HIImass[indxnew]/HIImass_exp, label = '%s'%(snap))

        if (len(indxnew) > 0):
                print('%s: min = %.4e, max = %.4e'%(snap, (HIImass[indxnew]/HIImass_exp).min(), (HIImass[indxnew]/HIImass_exp).max() ))

#ax.legend()

fig.savefig('HIIregion_test.png')

                                
