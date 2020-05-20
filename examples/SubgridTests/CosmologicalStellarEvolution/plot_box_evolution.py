###############################################################################
 # This file is part of SWIFT.
 # Copyright (c) 2018 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 # 
 # This program is free software: you can redistribute it and/or modify
 # it under the terms of the GNU Lesser General Public License as published
 # by the Free Software Foundation, either version 3 of the License, or
 # (at your option) any later version.
 # 
 # This program is distributed in the hope that it will be useful,
 # but WITHOUT ANY WARRANTY; without even the implied warranty of
 # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 # GNU General Public License for more details.
 # 
 # You should have received a copy of the GNU Lesser General Public License
 # along with this program.  If not, see <http://www.gnu.org/licenses/>.
 # 
 ##############################################################################

# Script used to plot time evolution of gas particle properties. Intended to 
# compare result of feedback due to one star placed in centre of uniform box
# of gas with output from EAGLE feedback test. Can also use as input output
# from SWIFT feedback test (tests/testFeedback) with the appropriate change
# to filepath.

import matplotlib
matplotlib.use("Agg")
from pylab import *
from scipy import stats
import h5py
import numpy as np
import glob
import os.path

# Plot parameters
params = {'axes.labelsize': 10,
'axes.titlesize': 10,
'font.size': 12,
'legend.fontsize': 12,
'xtick.labelsize': 10,
'ytick.labelsize': 10,
'text.usetex': True,
 'figure.figsize' : (12.90,6.45),
'figure.subplot.left'    : 0.07,
'figure.subplot.right'   : 0.995,
'figure.subplot.bottom'  : 0.06,
'figure.subplot.top'     : 0.92,
'figure.subplot.wspace'  : 0.3,
'figure.subplot.hspace'  : 0.2,
'lines.markersize' : 6,
'lines.linewidth' : 3.
}
rcParams.update(params)
rc('font',**{'family':'sans-serif','sans-serif':['Times']})


# Number of snapshots and elements
newest_snap_name = max(glob.glob('stellar_evolution_*.hdf5'))#, key=os.path.getctime)
n_snapshots = int(newest_snap_name.replace('stellar_evolution_','').replace('.hdf5','')) + 1
n_elements = 9

# Read the simulation data
sim = h5py.File("stellar_evolution_0000.hdf5", "r")
boxSize = sim["/Header"].attrs["BoxSize"][0]
time = sim["/Header"].attrs["Time"][0]
scheme = sim["/HydroScheme"].attrs["Scheme"]
kernel = sim["/HydroScheme"].attrs["Kernel function"]
neighbours = sim["/HydroScheme"].attrs["Kernel target N_ngb"]
eta = sim["/HydroScheme"].attrs["Kernel eta"]
git = sim["Code"].attrs["Git Revision"]
stellar_mass = sim["/PartType4/Masses"][0]
try: 
    E_SNII_cgs = double(sim["/Parameters"].attrs["EAGLEFeedback:SNII_energy_erg"])
    E_SNIa_cgs = double(sim["/Parameters"].attrs["EAGLEFeedback:SNIa_energy_erg"])
    ejecta_vel_cgs = double(sim["/Parameters"].attrs["EAGLEFeedback:AGB_ejecta_velocity_km_p_s"]) * 1e5
except:
    E_SNII_cgs = double(sim["/Parameters"].attrs["COLIBREFeedback:SNII_energy_erg"])
    E_SNIa_cgs = double(sim["/Parameters"].attrs["COLIBREFeedback:SNIa_energy_erg"])
    ejecta_vel_cgs = double(sim["/Parameters"].attrs["COLIBREFeedback:AGB_ejecta_velocity_km_p_s"]) * 1e5

try: 
    Eu_N_NSM_p_Msun = double(sim["/Parameters"].attrs["COLIBREFeedback:num_of_NSM_per_Msun"])
    Eu_N_CEJSN_p_Msun = double(sim["/Parameters"].attrs["COLIBREFeedback:num_of_CEJSN_per_Msun"])
    Eu_N_collapsar_p_Msun = double(sim["/Parameters"].attrs["COLIBREFeedback:num_of_collapsar_per_Msun"])
    Eu_yield_NSM_Msun = double(sim["/Parameters"].attrs["COLIBREFeedback:yield_Eu_from_NSM_event_Msun"])
    Eu_yield_CEJSN_Msun = double(sim["/Parameters"].attrs["COLIBREFeedback:yield_Eu_from_CEJSN_event_Msun"])
    Eu_yield_collapsar_Msun = double(sim["/Parameters"].attrs["COLIBREFeedback:yield_Eu_from_collapsar_event_Msun"])
except:
    Eu_N_NSM_p_Msun = 0.
    Eu_N_CEJSN_p_Msun = 0.
    Eu_N_collapsar_p_Msun = 0.
    Eu_yield_NSM_Msun = 0.
    Eu_yield_CEJSN_Msun = 0.
    Eu_yield_collapsar_Msun = 0.
    
# Units
unit_length_in_cgs = sim["/Units"].attrs["Unit length in cgs (U_L)"]
unit_mass_in_cgs = sim["/Units"].attrs["Unit mass in cgs (U_M)"]
unit_time_in_cgs = sim["/Units"].attrs["Unit time in cgs (U_t)"]
unit_temp_in_cgs = sim["/Units"].attrs["Unit temperature in cgs (U_T)"]
unit_vel_in_cgs = unit_length_in_cgs / unit_time_in_cgs
unit_energy_in_cgs = unit_mass_in_cgs * unit_vel_in_cgs * unit_vel_in_cgs
unit_density_in_cgs = unit_mass_in_cgs*unit_length_in_cgs**-3
unit_pressure_in_cgs = unit_mass_in_cgs/unit_length_in_cgs*unit_time_in_cgs**-2
unit_int_energy_in_cgs = unit_energy_in_cgs/unit_mass_in_cgs
unit_entropy_in_cgs = unit_energy_in_cgs/unit_temp_in_cgs
Gyr_in_cgs = sim["PhysicalConstants/CGS"].attrs["year"] * 1e9
Msun_in_cgs = sim["PhysicalConstants/CGS"].attrs["solar_mass"]

# Declare arrays to store SWIFT data
swift_box_gas_mass = zeros(n_snapshots)
swift_box_gas_mass_AGB = zeros(n_snapshots)
swift_box_gas_mass_SNII = zeros(n_snapshots)
swift_box_gas_mass_SNIa = zeros(n_snapshots)
swift_box_star_mass = zeros(n_snapshots)
swift_box_gas_metal_mass = zeros(n_snapshots)
swift_box_gas_metal_mass_AGB = zeros(n_snapshots)
swift_box_gas_metal_mass_SNII = zeros(n_snapshots)
swift_box_gas_metal_mass_SNIa = zeros(n_snapshots)
swift_element_mass = zeros((n_snapshots,n_elements))
swift_rprocess_mass = zeros(n_snapshots)
swift_box_gas_Eu_mass_NSM = zeros(n_snapshots)
swift_box_gas_Eu_mass_CEJSN = zeros(n_snapshots)
swift_box_gas_Eu_mass_Collapsar = zeros(n_snapshots)
swift_internal_energy = zeros(n_snapshots)
swift_kinetic_energy = zeros(n_snapshots)
swift_total_energy = zeros(n_snapshots)
swift_mean_u_start = 0.
t = zeros(n_snapshots)
z = zeros(n_snapshots)
a = zeros(n_snapshots)

# Read data from snapshots
for i in range(n_snapshots):
        #print("reading snapshot "+str(i))
        sim = h5py.File("stellar_evolution_%04d.hdf5"%i, "r")
        t[i] = sim["/Header"].attrs["Time"][0]
        z[i] = sim["/Cosmology"].attrs["Redshift"][0]
        a[i] = sim["/Cosmology"].attrs["Scale-factor"][0]
        
        masses = sim["/PartType0/Masses"][:]
        swift_box_gas_mass[i] = np.sum(masses)

        AGB_mass = sim["/PartType0/MassesFromAGB"][:]
        SNII_mass = sim["/PartType0/MassesFromSNII"][:]
        SNIa_mass = sim["/PartType0/MassesFromSNIa"][:]

        swift_box_gas_mass_AGB[i] = np.sum(AGB_mass)
        swift_box_gas_mass_SNII[i] = np.sum(SNII_mass)
        swift_box_gas_mass_SNIa[i] = np.sum(SNIa_mass)

        Z_star = sim["/PartType4/MetalMassFractions"][0]
        star_masses = sim["/PartType4/Masses"][:]
        swift_box_star_mass[i] = np.sum(star_masses)

        metallicities = sim["/PartType0/MetalMassFractions"][:]
        swift_box_gas_metal_mass[i] = np.sum(metallicities * masses)

        AGB_Z_frac = sim["/PartType0/MetalMassFractionsFromAGB"][:]
        SNII_Z_frac = sim["/PartType0/MetalMassFractionsFromSNII"][:]
        SNIa_Z_frac = sim["/PartType0/MetalMassFractionsFromSNIa"][:]

        swift_box_gas_metal_mass_AGB[i] = np.sum(AGB_Z_frac * masses)
        swift_box_gas_metal_mass_SNII[i] = np.sum(SNII_Z_frac * masses)
        swift_box_gas_metal_mass_SNIa[i] = np.sum(SNIa_Z_frac * masses)

        element_abundances = sim["/PartType0/ElementMassFractions"][:][:]
        for j in range(n_elements):
                swift_element_mass[i,j] = np.sum(element_abundances[:,j] * masses)

        try:
            swift_rprocess_mass[i] = np.sum(element_abundances[:,9] * masses)
            swift_box_gas_Eu_mass_NSM[i] = np.sum(sim["/PartType0/MassesFromNSM"][:])
            swift_box_gas_Eu_mass_CEJSN[i] = np.sum(sim["/PartType0/MassesFromCEJSN"][:])
            swift_box_gas_Eu_mass_Collapsar[i] = np.sum(sim["/PartType0/MassesFromCollapsar"][:])
        except:
            swift_rprocess_mass[i] = 0.
            swift_box_gas_Eu_mass_NSM[i] = 0.
            swift_box_gas_Eu_mass_CEJSN[i] = 0.
            swift_box_gas_Eu_mass_Collapsar[i] = 0.
            
        v = sim["/PartType0/Velocities"][:,:]
        v2 = v[:,0]**2 + v[:,1]**2 + v[:,2]**2
        u = sim["/PartType0/InternalEnergies"][:]
        swift_internal_energy[i] = np.sum(masses * u)
        swift_kinetic_energy[i] = np.sum(0.5 * masses * v2)
        swift_total_energy[i] = swift_kinetic_energy[i] + swift_internal_energy[i]

        if i == 0:
                swift_mean_u_start = np.mean(u)
        
        sim.close()

# Read expected yields from EAGLE. Choose which file to use based on metallicity used when
# running SWIFT (can be specified in yml file)
filename = "./StellarEvolutionSolution/Z_%.4f/StellarEvolutionTotal.txt"%Z_star

# Read EAGLE test output
data = loadtxt(filename)
eagle_time_Gyr = data[:,0]
eagle_total_mass = data[:,1] * stellar_mass / Msun_in_cgs * unit_mass_in_cgs
eagle_total_metal_mass = data[:,2] * stellar_mass / Msun_in_cgs * unit_mass_in_cgs
eagle_total_element_mass = data[:, 3:3+n_elements] * stellar_mass / Msun_in_cgs * unit_mass_in_cgs

eagle_energy_from_mass_cgs = eagle_total_mass * Msun_in_cgs * swift_mean_u_start * unit_int_energy_in_cgs
eagle_energy_ejecta_cgs = 0.5 * (eagle_total_mass * Msun_in_cgs) * ejecta_vel_cgs**2 

# Read the mass per channel
filename = "./StellarEvolutionSolution/Z_%.4f/StellarEvolutionAGB.txt"%Z_star
data = loadtxt(filename)
eagle_total_mass_AGB = data[:,1] * stellar_mass / Msun_in_cgs * unit_mass_in_cgs
eagle_total_metal_mass_AGB = data[:,2] * stellar_mass / Msun_in_cgs * unit_mass_in_cgs

filename = "./StellarEvolutionSolution/Z_%.4f/StellarEvolutionII.txt"%Z_star
data = loadtxt(filename)
eagle_total_mass_SNII = data[:,1] * stellar_mass / Msun_in_cgs * unit_mass_in_cgs
eagle_total_metal_mass_SNII = data[:,2] * stellar_mass / Msun_in_cgs * unit_mass_in_cgs

filename = "./StellarEvolutionSolution/Z_%.4f/StellarEvolutionIa.txt"%Z_star
data = loadtxt(filename)
eagle_total_mass_SNIa = data[:,1] * stellar_mass / Msun_in_cgs * unit_mass_in_cgs
eagle_total_metal_mass_SNIa = data[:,2] * stellar_mass / Msun_in_cgs * unit_mass_in_cgs

# Construct EAGLE scale-factors
t_Gyr = t * unit_time_in_cgs / Gyr_in_cgs
eagle_a = np.zeros(np.size(eagle_time_Gyr))
for i in range(np.size(eagle_a)):
    eagle_a[i] = np.interp(eagle_time_Gyr[i], t_Gyr, a)

# Construct Europium expectations
Eu_t = np.linspace(t[0], t[-1], 10000) * unit_time_in_cgs / Gyr_in_cgs
Eu_a = np.zeros(np.size(Eu_t))
for i in range(np.size(Eu_a)):
    Eu_a[i] = np.interp(Eu_t[i], t_Gyr, a)
Eu_mass_collapsar_Msun = zeros(np.size(Eu_t))
Eu_mass_CEJSN_Msun = zeros(np.size(Eu_t))
Eu_mass_NSM_Msun = zeros(np.size(Eu_t))

t0 = 0.03
t1 = 0.1
mask = Eu_t > t0
Eu_mass_NSM_Msun[mask] = Eu_N_NSM_p_Msun * (swift_box_star_mass[0] * unit_mass_in_cgs / Msun_in_cgs) * Eu_yield_NSM_Msun * np.log(Eu_t[mask] / t0)

mask = np.logical_and(Eu_t > t0, Eu_t <= t1)
Eu_mass_collapsar_Msun[mask] = Eu_N_collapsar_p_Msun * (swift_box_star_mass[0] * unit_mass_in_cgs / Msun_in_cgs)* Eu_yield_collapsar_Msun * (Eu_t[mask] - 0.03) / (0.1 - 0.03)
Eu_mass_CEJSN_Msun[mask] = Eu_N_CEJSN_p_Msun * (swift_box_star_mass[0] * unit_mass_in_cgs / Msun_in_cgs)* Eu_yield_CEJSN_Msun * (Eu_t[mask] - 0.03) / (0.1 - 0.03)

mask = Eu_t > t1
Eu_mass_collapsar_Msun[mask] = Eu_N_collapsar_p_Msun * (swift_box_star_mass[0] * unit_mass_in_cgs / Msun_in_cgs)* Eu_yield_collapsar_Msun
Eu_mass_CEJSN_Msun[mask] = Eu_N_CEJSN_p_Msun * (swift_box_star_mass[0] * unit_mass_in_cgs / Msun_in_cgs)* Eu_yield_CEJSN_Msun

Eu_mass_total_Msun = Eu_mass_collapsar_Msun + Eu_mass_CEJSN_Msun + Eu_mass_NSM_Msun

# Plot the interesting quantities
figure()

redshift_ticks = np.array([0., 0.2, 0.5, 1., 2., 3., 5., 10.])
redshift_labels = ["$0$", "$0.2$", "$0.5$", "$1$", "$2$", "$3$", "$5$", "$10$"]
a_ticks = 1. / (redshift_ticks + 1.)

suptitle("Star metallicity Z = %.4f"%Z_star)

# Box gas mass --------------------------------
subplot(231, xscale="log")
plot(a[1:], (swift_box_gas_mass[1:] - swift_box_gas_mass[0])* unit_mass_in_cgs / Msun_in_cgs, linewidth=0.5, color='k', label='Total')
plot(a, swift_box_gas_mass_AGB * unit_mass_in_cgs / Msun_in_cgs, linewidth=0.5, color='C0', label='AGB')
plot(a, swift_box_gas_mass_SNII * unit_mass_in_cgs / Msun_in_cgs, linewidth=0.5, color='C1', label='SNII')
plot(a, swift_box_gas_mass_SNIa * unit_mass_in_cgs / Msun_in_cgs, linewidth=0.5, color='C2', label='SNIa')
plot(eagle_a[1:],eagle_total_mass[:-1],linewidth=0.5, color='k', ls='--')
plot(eagle_a[1:],eagle_total_mass_AGB[:-1],linewidth=0.5, color='C0', ls='--')
plot(eagle_a[1:],eagle_total_mass_SNII[:-1],linewidth=0.5, color='C1', ls='--')
plot(eagle_a[1:],eagle_total_mass_SNIa[:-1],linewidth=0.5, color='C2', ls='--')
legend(loc="upper left", ncol=2, fontsize=8)
xlabel("${\\rm Redshift}$", labelpad=0)
ylabel("Change in total gas particle mass ${[\\rm M_\\odot]}$", labelpad=2)
xlim(0.08, 1.05)
ticklabel_format(style='sci', axis='y', scilimits=(0,0))
xticks(a_ticks, redshift_labels)
gca().tick_params(axis='x', which='minor', bottom=False)


# Box star mass --------------------------------
subplot(232, xscale="log")
plot(a, (swift_box_star_mass)* unit_mass_in_cgs / Msun_in_cgs, linewidth=0.5, color='k', label='SWIFT')
plot(eagle_a[1:], swift_box_star_mass[0] * unit_mass_in_cgs / Msun_in_cgs - eagle_total_mass[:-1],linewidth=0.5,color='k',label='EAGLE test', ls='--')
xlabel("${\\rm Redshift}$", labelpad=0)
ylabel("Change in total star particle mass ${[\\rm M_\\odot]}$", labelpad=2)
ticklabel_format(style='sci', axis='y', scilimits=(0,0))
xlim(0.08, 1.05)
xticks(a_ticks, redshift_labels)
gca().tick_params(axis='x', which='minor', bottom=False)
legend()

# Box gas element  mass --------------------------------
colours = ['k','r','g','b','c','y','m','skyblue','plum']
element_names = ['H','He','C','N','O','Ne','Mg','Si','Fe']
subplot(233, xscale="log", yscale="log")
for j in range(n_elements):
        plot(a[1:], (swift_element_mass[1:,j] - swift_element_mass[0,j]) * unit_mass_in_cgs / Msun_in_cgs, linewidth=0.5, color=colours[j], ms=0.5, label=element_names[j])
        plot(eagle_a[1:],eagle_total_element_mass[:-1,j],linewidth=1,color=colours[j],linestyle='--')
xlabel("${\\rm Redshift}$", labelpad=0)
ylabel("Change in element mass of gas particles ${[\\rm M_\\odot]}$", labelpad=2)
legend(bbox_to_anchor=(-0.135, 1.), ncol=1, fontsize=8, handlelength=1)
xlim(0.08, 1.05)
ylim(3e-4 * swift_box_star_mass[0] * unit_mass_in_cgs / Msun_in_cgs, swift_box_star_mass[0] * unit_mass_in_cgs / Msun_in_cgs)
xticks(a_ticks, redshift_labels)
gca().tick_params(axis='x', which='minor', bottom=False)

subplot(234, xscale="log")
plot(a, swift_rprocess_mass* unit_mass_in_cgs / Msun_in_cgs, linewidth=0.5, color='k', label="Total")
plot(a, swift_box_gas_Eu_mass_NSM * unit_mass_in_cgs / Msun_in_cgs, linewidth=0.5, color='C2', label='NSM')
plot(a, swift_box_gas_Eu_mass_CEJSN * unit_mass_in_cgs / Msun_in_cgs, linewidth=0.5, color='C1', label='CEJSN')
plot(a, swift_box_gas_Eu_mass_Collapsar * unit_mass_in_cgs / Msun_in_cgs, linewidth=0.5, color='C0', label='Collapsars')
plot(Eu_a, Eu_mass_total_Msun, linewidth=0.5, color='k', ls='--')
plot(Eu_a, Eu_mass_collapsar_Msun, linewidth=0.5, color='C0', ls='--')
plot(Eu_a, Eu_mass_CEJSN_Msun, linewidth=0.5, color='C1', ls='--')
plot(Eu_a, Eu_mass_NSM_Msun, linewidth=0.5, color='C2', ls='--')
legend(loc="upper left", ncol=2, fontsize=8)
xlabel("${\\rm Redshift}$", labelpad=0)
ylabel("Change in Europium mass of gas particles ${[\\rm M_\\odot]}$", labelpad=2)
xlim(0.08, 1.05)
xticks(a_ticks, redshift_labels)
gca().tick_params(axis='x', which='minor', bottom=False)

# Box gas metal mass --------------------------------
subplot(235, xscale="log")
plot(a[1:] , (swift_box_gas_metal_mass[1:] - swift_box_gas_metal_mass[0])* unit_mass_in_cgs / Msun_in_cgs, linewidth=0.5, color='k', label='Total')
plot(a, swift_box_gas_metal_mass_AGB * unit_mass_in_cgs / Msun_in_cgs, linewidth=0.5, color='C0', label='AGB')
plot(a, swift_box_gas_metal_mass_SNII * unit_mass_in_cgs / Msun_in_cgs, linewidth=0.5, color='C1', label='SNII')
plot(a, swift_box_gas_metal_mass_SNIa * unit_mass_in_cgs / Msun_in_cgs, linewidth=0.5, color='C2', label='SNIa')
plot(eagle_a[1:],eagle_total_metal_mass[:-1],linewidth=0.5,color='k', ls='--')
plot(eagle_a[1:],eagle_total_metal_mass_AGB[:-1],linewidth=0.5, color='C0', ls='--')
plot(eagle_a[1:],eagle_total_metal_mass_SNII[:-1],linewidth=0.5, color='C1', ls='--')
plot(eagle_a[1:],eagle_total_metal_mass_SNIa[:-1],linewidth=0.5, color='C2', ls='--')
legend(loc="center right", ncol=2, fontsize=8)
xlabel("${\\rm Redshift}$", labelpad=0)
ylabel("Change in total metal mass of gas particles ${[\\rm M_\\odot]}$", labelpad=2)
xlim(0.08, 1.05)
xticks(a_ticks, redshift_labels)
gca().tick_params(axis='x', which='minor', bottom=False)
ticklabel_format(style='sci', axis='y', scilimits=(0,0))

savefig("box_evolution_Z_%.4f.png"%(Z_star), dpi=200)

