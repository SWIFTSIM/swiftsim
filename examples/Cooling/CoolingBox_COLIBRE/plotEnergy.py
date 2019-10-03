from h5py import File
import numpy as np
import matplotlib
from glob import glob
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Plot parameters
params = {
    'axes.labelsize': 10,
    'axes.titlesize': 10,
    'font.size': 12,
    'legend.fontsize': 12,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'text.usetex': True,
    'figure.figsize': (5, 5),
    'figure.subplot.left': 0.1,
    'figure.subplot.right': 0.99,
    'figure.subplot.bottom': 0.1,
    'figure.subplot.top': 0.96,
    'figure.subplot.wspace': 0.15,
    'figure.subplot.hspace': 0.12,
    'lines.markersize': 6,
    'lines.linewidth': 3.,
}
plt.rcParams.update(params)


# Some constants in cgs units
k_b_cgs = 1.38064852e-16   # boltzmann
m_h_cgs = 1.672621898e-24  # proton mass


# File containing the total energy
stats_filename = "./energy.txt"

# First snapshot
snap_filename = "output_0000.hdf5"

# Read the initial state of the gas
f = File(snap_filename, 'r')

# Read the units parameters from the snapshot
units = f["InternalCodeUnits"]
UnitMass_in_cgs = units.attrs["Unit mass in cgs (U_M)"]
UnitLength_in_cgs = units.attrs["Unit length in cgs (U_L)"]
UnitTime_in_cgs = units.attrs["Unit time in cgs (U_t)"]

# Read the adiabatic index
gamma = float(f["HydroScheme"].attrs["Adiabatic index"])

# Unit conversions 
convert_seconds_to_Myr = 1.0 / (60. * 60. * 24. * 365.25 * 1.0e6) 

def convert_energy_to_cgs(u):
    """ Convert total energy from code units to erg."""
    return u * UnitMass_in_cgs * ((UnitLength_in_cgs / UnitTime_in_cgs) ** 2.0) 

# Read energy and time arrays
array = np.genfromtxt(stats_filename, skip_header=1)
time = array[:, 0] * UnitTime_in_cgs * convert_seconds_to_Myr
total_mass = array[:, 1]
total_energy = array[:, 2]
kinetic_energy = array[:, 3]
internal_energy = array[:, 4]
radiated_energy = array[:, 8]
initial_energy = total_energy[0]

# Conversions to cgs
total_energy_cgs = convert_energy_to_cgs(total_energy) 
kinetic_energy_cgs = convert_energy_to_cgs(kinetic_energy)
internal_energy_cgs = convert_energy_to_cgs(internal_energy)
radiated_energy_cgs = convert_energy_to_cgs(radiated_energy)

# Read snapshots
files = glob("output_*.hdf5")
N = len(files)
u_snap_cgs = np.zeros(N)
time_snap_cgs = np.zeros(N)
for i in range(N):
    snap = File(files[i], 'r')
    u = snap["/PartType0/InternalEnergies"][:] * snap["/PartType0/Masses"][:]
    u = sum(u) 
    u_snap_cgs[i] = convert_energy_to_cgs(u)
    time_snap_cgs[i] = snap["/Header"].attrs["Time"] * UnitTime_in_cgs * convert_seconds_to_Myr


plt.figure()

# Total energy from snapshots 
plt.plot(time_snap_cgs, u_snap_cgs, 'rD', markersize = 3)

# Energies from Stats file 
plt.plot(time, total_energy_cgs, 'r-', linewidth = 1.2, label="Gas total energy")
plt.plot(time, radiated_energy_cgs, 'g-', linewidth = 1.2, label="Radiated energy")
plt.plot(time, total_energy_cgs + radiated_energy_cgs, 'b-', linewidth = 1.2, label="Gas total + radiated")

plt.legend(loc="right", fontsize=8, frameon=False,
           handlelength=3, ncol=1)
plt.xlabel("${\\rm{Time~[Myr]}}$", labelpad = 0)
plt.ylabel("${\\rm{Energy~[erg]}}$", labelpad = 0)

plt.savefig("energy.png", dpi=200)
