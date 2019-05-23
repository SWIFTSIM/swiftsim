#!/usr/bin/env python3
################################################################################
# This file is part of SWIFT.
# Copyright (c) 2019 Folkert Nobels (nobels@strw.leidenuniv.nl)
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
################################################################################
import matplotlib
matplotlib.use("Agg")
import numpy as np
import swiftsimio as sw
import matplotlib.pyplot as plt
import unyt
import h5py as h5
import sys
from scipy import stats
from sphviewer.tools import QuickView
from matplotlib.colors import LogNorm

# Plot parameters
params = {
    "axes.labelsize": 10,
    "axes.titlesize": 10,
    "font.size": 9,
    "legend.fontsize": 9,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "text.usetex": True,
    "figure.figsize": (3.15, 3.15),
    "figure.subplot.left": 0.15,
    "figure.subplot.right": 0.99,
    "figure.subplot.bottom": 0.13,
    "figure.subplot.top": 0.99,
    "figure.subplot.wspace": 0.15,
    "figure.subplot.hspace": 0.12,
    "lines.markersize": 6,
    "lines.linewidth": 2.0,
    "text.latex.unicode": True,
}
plt.rcParams.update(params)
plt.rc("font", **{"family": "sans-serif", "sans-serif": ["Times"]})

snap = int(sys.argv[1])
filename = "output_%.4d.hdf5"%snap

f = h5.File(filename, "r")


# Physical constants
k_in_cgs = 1.38064852e-16
mH_in_cgs = 1.6737236e-24
year_in_cgs = 3600.0 * 24 * 365.0
Msun_in_cgs = 1.98848e33
G_in_cgs = 6.67259e-8
pc_in_cgs = 3.08567758e18
Msun_p_pc2 = Msun_in_cgs / pc_in_cgs**2

# Read units 
unit_length_in_cgs = f["/Units"].attrs["Unit length in cgs (U_L)"] 
unit_mass_in_cgs = f["/Units"].attrs["Unit mass in cgs (U_M)"] 
unit_time_in_cgs = f["/Units"].attrs["Unit time in cgs (U_t)"]

boxsize = f["/Header"].attrs["BoxSize"]
centre = boxsize / 2.0

# Read gas properties
gas_pos = f["/PartType0/Coordinates"][:, :]
gas_mass = f["/PartType0/Masses"][:]
gas_rho = f["/PartType0/Density"][:]
gas_T = f["/PartType0/Temperature"][:]
gas_SFR = f["/PartType0/SFR"][:]
gas_XH = f["/PartType0/ElementAbundance"][:, 0]
gas_Z = f["/PartType0/Metallicity"][:]
gas_hsml = f["/PartType0/SmoothingLength"][:]

# Read the Star properties
stars_pos = f["/PartType4/Coordinates"][:, :]
stars_BirthDensity = f["/PartType4/BirthDensity"][:]
stars_BirthTime = f["/PartType4/BirthTime"][:]
stars_XH = f["/PartType4/ElementAbundance"][:,0]

# Centre the box
gas_pos[:, 0] -= centre[0]
gas_pos[:, 1] -= centre[1]
gas_pos[:, 2] -= centre[2]

stars_pos[:,0] -= centre[0]
stars_pos[:,1] -= centre[1]
stars_pos[:,2] -= centre[2]

# Calculate Gravitational constant in internal units
G = G_in_cgs * ( unit_length_in_cgs**3 / unit_mass_in_cgs / unit_time_in_cgs**2)**(-1)

# Read the paramters of the SF model
f_star = float(f["/Parameters"].attrs["COLIBREStarFormation:star_formation_efficiency"])
temp_thresh = float(f["/Parameters"].attrs["COLIBREStarFormation:temperature_threshold"])
max_den = float(f["/Parameters"].attrs["COLIBREStarFormation:threshold_max_density_H_p_cm3"])
T_dex = float(f["/Parameters"].attrs["COLIBREStarFormation:KS_temperature_margin_dex"])

# Read parameters of the entropy floor
EAGLEfloor_Jeans_rho_norm = float(f["/Parameters"].attrs["EAGLEEntropyFloor:Jeans_density_threshold_H_p_cm3"])
EAGLEfloor_Jeans_temperature_norm_K = float(f["/Parameters"].attrs["EAGLEEntropyFloor:Jeans_temperature_norm_K"])
EAGLEfloor_Jeans_gamma_effective = float(f["/Parameters"].attrs["EAGLEEntropyFloor:Jeans_gamma_effective"])
EAGLEfloor_cool_rho_norm = float(f["/Parameters"].attrs["EAGLEEntropyFloor:Cool_density_threshold_H_p_cm3"])
EAGLEfloor_cool_temperature_norm_K = float(f["/Parameters"].attrs["EAGLEEntropyFloor:Cool_temperature_norm_K"])
EAGLEfloor_cool_gamma_effective = float(f["/Parameters"].attrs["EAGLEEntropyFloor:Cool_gamma_effective"])

# Read reference metallicity
#EAGLE_Z = float(f["/Parameters"].attrs["EAGLEChemistry:init_abundance_metal"])

# Turn the mass into better units
gas_mass *= unit_mass_in_cgs / Msun_in_cgs

# Turn the SFR into better units
gas_SFR = np.maximum(gas_SFR, np.zeros(np.size(gas_SFR)))
gas_SFR /= unit_time_in_cgs / year_in_cgs
gas_SFR *= unit_mass_in_cgs / Msun_in_cgs

# Make it a Hydrogen number density
gas_nH = gas_rho * unit_mass_in_cgs / unit_length_in_cgs ** 3
gas_nH /= mH_in_cgs
gas_nH *= gas_XH

stars_BirthDensity *= unit_mass_in_cgs / unit_length_in_cgs ** 3
stars_BirthDensity /= mH_in_cgs
stars_BirthDensity *= stars_XH

gas_sSFR = gas_SFR / gas_mass

# Calculate the median gas mass
median_gas_mass = np.median(gas_mass)

# Equations of state
eos_cool_rho = np.logspace(-5, 5, 1000)
eos_cool_T = EAGLEfloor_cool_temperature_norm_K * (eos_cool_rho / EAGLEfloor_cool_rho_norm) ** ( EAGLEfloor_cool_gamma_effective - 1.0 )
eos_Jeans_rho = np.logspace(-1, 5, 1000)
eos_Jeans_T = EAGLEfloor_Jeans_temperature_norm_K * (eos_Jeans_rho / EAGLEfloor_Jeans_rho_norm) ** (EAGLEfloor_Jeans_gamma_effective - 1.0 ) 

#################################################################

# Plot the phase space diagram
plt.figure()
plt.subplot(111, xscale="log", yscale="log")
plt.plot(eos_cool_rho, eos_cool_T, "k--", lw=0.6)
plt.plot(eos_Jeans_rho, eos_Jeans_T, "k--", lw=0.6)
plt.scatter(gas_nH, gas_T, s=0.2)
plt.xlabel("${\\rm Density}~n_{\\rm H}~[{\\rm cm^{-3}}]$", labelpad=0)
plt.ylabel("${\\rm Temperature}~T~[{\\rm K}]$", labelpad=2)
plt.xlim(1e-4, 3e6)
plt.ylim(10.0, 2e5)
plt.savefig("rhoT.png", dpi=200)
plt.close()

# Plot the phase space diagram for SF gas
plt.figure()
plt.subplot(111, xscale="log", yscale="log")
plt.plot(eos_cool_rho, eos_cool_T, "k--", lw=0.6)
plt.plot(eos_Jeans_rho, eos_Jeans_T, "k--", lw=0.6)
plt.scatter(gas_nH[gas_SFR > 0.0], gas_T[gas_SFR > 0.0], s=0.2)
# Make the temperature threshold line
max_eos_cool_rho = np.min(eos_Jeans_rho)
if T_dex<10.:
    jeans_T = eos_Jeans_T*10**T_dex
    plt.plot(eos_cool_rho[eos_cool_rho<max_eos_cool_rho], eos_cool_T[eos_cool_rho<max_eos_cool_rho]*10**T_dex, "k--", lw=0.5)
    plt.plot(eos_Jeans_rho, jeans_T, "--",color='gray', lw=0.4)
    plt.plot(eos_Jeans_rho[jeans_T < temp_thresh], jeans_T[jeans_T < temp_thresh], "k--", lw=0.5)
    plt.axhline(y=temp_thresh,linestyle="--",color='gray', lw=0.4)
    plt.plot([eos_Jeans_rho[jeans_T < temp_thresh][-1],3e6],[temp_thresh,temp_thresh], "k--", lw=0.5)
    l2 = np.array((eos_cool_rho[eos_cool_rho<max_eos_cool_rho][-1],eos_cool_T[eos_cool_rho<max_eos_cool_rho][-1]*10**T_dex))
    trans_angle = plt.gca().transData.transform_angles(np.array((30,)),
                                                   l2.reshape((1, 2)))
else: 
    plt.axhline(y=temp_thresh,linestyle="--",color='k', lw=0.4)

plt.text(1e2,temp_thresh,"Temperature threshold", fontsize=8, ha="right", va="bottom")
#plt.text(l2[0],l2[1],"Temperature threshold!", rotation=trans_angle,rotation_mode="anchor")
plt.xlabel("${\\rm Density}~n_{\\rm H}~[{\\rm cm^{-3}}]$", labelpad=0)
plt.ylabel("${\\rm Temperature}~T~[{\\rm K}]$", labelpad=2)
plt.xlim(1e-4, 3e6)
plt.ylim(10.0, 2e5)
plt.savefig("rhoT_SF.png", dpi=200)
plt.close()

def sSFR_COLIBRE(rho,fstar):
    return fstar*(32*rho*G_in_cgs/(3*np.pi))**.5

rho = np.logspace(1,7,1000)
sSFR_theory = sSFR_COLIBRE(rho*1.6737236e-24*1.3,f_star)
sSFR_theory *= 365*24*3600 * 1e9   

# density vs sSFR
plt.rcParams.update({"figure.subplot.left": 0.18})
plt.figure()
plt.subplot(111, xscale="log", yscale="log")
plt.scatter(gas_nH, gas_sSFR*1e9,s=1,label='Simulation')
plt.xlabel("Number density [$\\rm cm^{-3}$]")
plt.ylabel('Specific star formation rate [$\\rm Gyr^{-1}$]')
plt.plot(rho,sSFR_theory,'--r',label='sSFR COLIBRE',linewidth=0.5)
plt.savefig("density-sSFR.png",dpi=200)
plt.close()
plt.rcParams.update({"figure.subplot.left": 0.15})

SFR_theory = 10**(np.log10(sSFR_theory)+np.log10(median_gas_mass)-9.)

# 3D Density vs SFR
plt.rcParams.update({"figure.subplot.left": 0.18})
plt.figure()
plt.subplot(111, xscale="log", yscale="log")
plt.scatter(gas_nH, gas_SFR, s=0.2)
plt.plot(rho,SFR_theory,'k--',lw=1,label='SFR COLIBRE')
plt.xlabel("${\\rm Density}~n_{\\rm H}~[{\\rm cm^{-3}}]$", labelpad=0)
plt.ylabel("${\\rm SFR}~[{\\rm M_\\odot~\\cdot~yr^{-1}}]$", labelpad=2)
#plt.xlim(1e-2, 1e5)
#plt.ylim(10**SFR_low_min, 10**(SFR_high_max+0.1))
plt.savefig("rho_SFR.png", dpi=200)
plt.close()
plt.rcParams.update({"figure.subplot.left": 0.15})

########################################################################3

star_mask = (
    (stars_pos[:, 0] > -15)
    & (stars_pos[:, 0] < 15)
    & (stars_pos[:, 1] > -15)
    & (stars_pos[:, 1] < 15)
    & (stars_pos[:, 2] < 1.0)
    & (stars_pos[:, 2] > -1.0)
)

stars_BirthDensity = stars_BirthDensity[star_mask] 
#stars_BirthFlag = stars_BirthFlag[star_mask]
stars_BirthTime = stars_BirthTime[star_mask]

# Histogram of the birth density
plt.figure()
plt.subplot(111, xscale="linear", yscale="linear")
plt.hist(np.log10(stars_BirthDensity),density=True,bins=20,range=[-2,5])
plt.xlabel("${\\rm Stellar~birth~density}~n_{\\rm H}~[{\\rm cm^{-3}}]$", labelpad=0)
plt.ylabel("${\\rm Probability}$", labelpad=-7)
plt.savefig("BirthDensity.png", dpi=200)
plt.close()

########################################################################3

# Select gas in a pillow box around the galaxy
mask = (
    (gas_pos[:, 0] > -15)
    & (gas_pos[:, 0] < 15)
    & (gas_pos[:, 1] > -15)
    & (gas_pos[:, 1] < 15)
    & (gas_pos[:, 2] < 1.0)
    & (gas_pos[:, 2] > -1.0)
)
gas_pos = gas_pos[mask, :]
gas_SFR = gas_SFR[mask]
gas_nH = gas_nH[mask]
gas_rho = gas_rho[mask]
gas_T = gas_T[mask]
gas_mass = gas_mass[mask]
gas_Z = gas_Z[mask]
gas_hsml = gas_hsml[mask]


# Make a crude map of the gas
plt.figure()
plt.subplot(111)
plt.scatter(gas_pos[:, 0], gas_pos[:, 1], s=0.1)
plt.xlabel("${\\rm Pos}~x~[{\\rm kpc}]$", labelpad=0)
plt.ylabel("${\\rm Pos}~y~[{\\rm kpc}]$", labelpad=-3)
plt.xlim(-12, 12)
plt.ylim(-12, 12)
plt.savefig("face_on.png", dpi=200)
plt.close()

plt.figure()
plt.subplot(111)
plt.scatter(gas_pos[:, 0], gas_pos[:, 2], s=0.1)
plt.xlabel("${\\rm Pos}~x~[{\\rm kpc}]$", labelpad=0)
plt.ylabel("${\\rm Pos}~z~[{\\rm kpc}]$", labelpad=-3)
plt.xlim(-12, 12)
plt.ylim(-12, 12)
plt.savefig("edge_on.png", dpi=200)
plt.close()

# Now a SF map
plt.rcParams.update({"figure.figsize": (4.15, 3.15)})
plt.figure()
plt.subplot(111)
plt.scatter(gas_pos[:, 0], gas_pos[:, 1], s=0.1, c=gas_SFR)
plt.xlabel("${\\rm Pos}~x~[{\\rm kpc}]$", labelpad=0)
plt.ylabel("${\\rm Pos}~y~[{\\rm kpc}]$", labelpad=-3)
plt.colorbar()
plt.xlim(-12, 12)
plt.ylim(-12, 12)
plt.savefig("SF_face_on.png", dpi=200)
plt.close()

########################################################################3

# Bin the data in kpc-size patches

x_edges = np.linspace(-15, 15, 31)
y_edges = np.linspace(-15, 15, 31)

map_mass, _, _, _ = stats.binned_statistic_2d(
    gas_pos[:, 0], gas_pos[:, 1], gas_mass, statistic="sum", bins=(x_edges, y_edges)
)
map_SFR, _, _, _ = stats.binned_statistic_2d(
    gas_pos[:, 0], gas_pos[:, 1], gas_SFR, statistic="sum", bins=(x_edges, y_edges)
)

# Mass map
plt.figure()
plt.subplot(111)
plt.pcolormesh(x_edges, y_edges, np.log10(map_mass))
plt.colorbar()
plt.xlim(-12, 12)
plt.ylim(-12, 12)
plt.xlabel("${\\rm Pos}~x~[{\\rm kpc}]$", labelpad=0)
plt.ylabel("${\\rm Pos}~y~[{\\rm kpc}]$", labelpad=-3)
plt.savefig("Map_mass.png", dpi=200)
plt.close()

# SF map
plt.figure()
plt.subplot(111)
plt.pcolormesh(x_edges, y_edges, np.log10(map_SFR), vmax=-0.5, vmin=-4.5)
plt.colorbar()
plt.xlim(-12, 12)
plt.ylim(-12, 12)
plt.xlabel("${\\rm Pos}~x~[{\\rm kpc}]$", labelpad=0)
plt.ylabel("${\\rm Pos}~y~[{\\rm kpc}]$", labelpad=-3)
plt.savefig("Map_SFR.png", dpi=200)
plt.close()

#########################################################################

# Give a minimum SF surface density for the plots
map_SFR[map_SFR < 1e-6] = 1e-6

# Emperical KS law
def KS(sigma_g,n,A):
    return A*sigma_g**n

Sigma_g = np.logspace(1,3,1000)

Sigma_star = KS(Sigma_g,1.4,1.515e-4)

# KS relation
plt.rcParams.update({"figure.figsize": (3.15, 3.15), "figure.subplot.left": 0.18})
plt.figure()
plt.subplot(111, xscale="log", yscale="log")
plt.scatter(map_mass.flatten() / 1e6, map_SFR.flatten(), s=0.4)
plt.plot(Sigma_g,Sigma_star,'--r',linewidth=.5)
plt.xlim(0.3, 900)
plt.ylim(3e-7, 3e1)
plt.xlabel("$\\Sigma_g~[{\\rm M_\\odot\\cdot pc^{-2}}]$", labelpad=0)
plt.ylabel(
    "$\\Sigma_{\\rm SFR}~[{\\rm M_\\odot \\cdot yr^{-1} \\cdot kpc^{-2}}]$", labelpad=0
)
plt.savefig("KS_law.png", dpi=200)
plt.close()
plt.rcParams.update({"figure.figsize": (3.15, 3.15), "figure.subplot.left": 0.15})

try:
    # Read the logger file
    logdata = np.loadtxt("SFR.txt")

    # Define the logger data in the correct units
    timelog = logdata[:, 1] * 9.778131e2

    plt.plot(timelog, logdata[:, 7] * 1.023009e01, label="SFH log file")
    plt.xlabel("Time (Myr)")
    plt.ylabel("SFH [$\\rm M_\odot \\rm yr^{-1}$]")
    #plt.ylim(0, 15)
    #plt.xlim(0, 100)
    #plt.legend()
    plt.savefig("SFH_log_file.png", dpi=200)
    plt.close()
except: 
    print("SFH file not properly formatted!! Did you crash?")

qv = QuickView(
    gas_pos,
    mass=gas_mass,
    r="infinity",
    xsize=len(x_edges)-1,
    ysize=len(y_edges)-1,
    p=0,  # Viewing angle theta
    roll=0,  # Viewing angle phi
    plot=False,
    logscale=False,
    hsml=2*gas_hsml
)

map_mass2 = qv.get_image()
extent_mass = qv.get_extent()

gas_SFR[gas_SFR<=0] = 1e-10

qv = QuickView(
    gas_pos,
    mass=gas_SFR,
    r="infinity",
    xsize=len(x_edges)-1,
    ysize=len(y_edges)-1,
    p=0,  # Viewing angle theta
    roll=0,  # Viewing angle phi
    plot=False,
    logscale=False,
    hsml=2*gas_hsml
)

map_SFR2 = qv.get_image()
extent_SFR = qv.get_extent()

# Mass map 2
plt.figure()
plt.subplot(111)
plt.imshow(np.log10(map_mass2),extent=extent_mass)
plt.colorbar()
plt.xlim(-12, 12)
plt.ylim(-12, 12)
plt.xlabel("${\\rm Pos}~x~[{\\rm kpc}]$", labelpad=0)
plt.ylabel("${\\rm Pos}~y~[{\\rm kpc}]$", labelpad=-3)
plt.savefig("Map_mass_SPHVIEWER.png", dpi=200)

plot_map_SFR2 = np.zeros(np.shape(map_SFR2))
plot_map_SFR2[map_SFR2>1e-8] = map_SFR2[map_SFR2>1e-8]
# SF map 2
plt.figure()
plt.subplot(111)
plt.imshow(np.log10(plot_map_SFR2),extent=extent_SFR, vmax = -.5, vmin=-4.5)
plt.colorbar()
plt.xlim(-12, 12)
plt.ylim(-12, 12)
plt.xlabel("${\\rm Pos}~x~[{\\rm kpc}]$", labelpad=0)
plt.ylabel("${\\rm Pos}~y~[{\\rm kpc}]$", labelpad=-3)
plt.savefig("Map_SFR_SPHVIEWER.png", dpi=200)

# Make the KS plot using SPHviewer
plot_map_SFR2[plot_map_SFR2<=0] = 1e-10

plt.rcParams.update({    "figure.figsize": (3.15, 3.15),     "figure.subplot.left": 0.18,})
plt.figure()
plt.subplot(111, xscale="log", yscale="log")
plt.scatter(map_mass2.flatten() / 1e6, plot_map_SFR2.flatten(), s=0.4)
plt.plot(Sigma_g,Sigma_star,'--r',linewidth=.5)
plt.xlim(0.3, 900)
plt.ylim(3e-7, 3)
plt.xlabel("$\\Sigma_g~[{\\rm M_\\odot\\cdot pc^{-2}}]$", labelpad=0)
plt.ylabel("$\\Sigma_{\\rm SFR}~[{\\rm M_\\odot \\cdot yr^{-1} \\cdot kpc^{-2}}]$", labelpad=0)
plt.savefig("KS_law_SPHVIEWER.png", dpi=200)

############################################################################
# Bin the data in kpc-size patches



x_edges = np.linspace(-15, 15, 301)
y_edges = np.linspace(-15, 15, 301)

map_mass, _, _, _ = stats.binned_statistic_2d(
    gas_pos[:, 0], gas_pos[:, 1], gas_mass, statistic="sum", bins=(x_edges, y_edges)
)
map_SFR, _, _, _ = stats.binned_statistic_2d(
    gas_pos[:, 0], gas_pos[:, 1], gas_SFR, statistic="sum", bins=(x_edges, y_edges)
)

# Mass map
plt.rcParams.update({"figure.figsize": (3.15, 3.15), "figure.subplot.left": 0.18})
plt.figure()
plt.subplot(111)
plt.pcolormesh(x_edges, y_edges, np.log10(map_mass))
plt.colorbar()
plt.xlim(-12, 12)
plt.ylim(-12, 12)
plt.xlabel("${\\rm Pos}~x~[{\\rm kpc}]$", labelpad=0)
plt.ylabel("${\\rm Pos}~y~[{\\rm kpc}]$", labelpad=-3)
plt.savefig("Map_mass_subkpc.png", dpi=200)
plt.close()

# SF map
plt.figure()
plt.subplot(111)
plt.pcolormesh(x_edges, y_edges, np.log10(map_SFR)) #, vmax=-0.5, vmin=-4.5)
plt.colorbar()
plt.xlim(-12, 12)
plt.ylim(-12, 12)
plt.xlabel("${\\rm Pos}~x~[{\\rm kpc}]$", labelpad=0)
plt.ylabel("${\\rm Pos}~y~[{\\rm kpc}]$", labelpad=-3)
plt.savefig("Map_SFR_subkpc.png", dpi=200)
plt.close()

#########################################################################

# Give a minimum SF surface density for the plots
map_SFR *= 1e2 # correct for 10 times smaller scale.
map_SFR[map_SFR < 1e-6] = 1e-6

# Emperical KS law
def KS(sigma_g,n,A):
    return A*sigma_g**n

Sigma_g = np.logspace(1,3,1000)

Sigma_star = KS(Sigma_g,1.4,1.515e-4)

# KS relation
plt.rcParams.update({"figure.figsize": (3.15, 3.15), "figure.subplot.left": 0.15})
plt.figure()
plt.subplot(111, xscale="log", yscale="log")
plt.scatter(map_mass.flatten() / 1e6*1e2, map_SFR.flatten(), s=0.4)
plt.plot(Sigma_g,Sigma_star,'--r',linewidth=.5)
plt.xlim(0.3, 900)
plt.ylim(3e-7, 3e1)
plt.xlabel("$\\Sigma_g~[{\\rm M_\\odot\\cdot pc^{-2}}]$", labelpad=0)
plt.ylabel(
    "$\\Sigma_{\\rm SFR}~[{\\rm M_\\odot \\cdot yr^{-1} \\cdot kpc^{-2}}]$", labelpad=0
)
plt.savefig("KS_law_subkpc.png", dpi=200)
plt.close()
plt.rcParams.update({"figure.figsize": (3.15, 3.15), "figure.subplot.left": 0.18})

qv = QuickView(
    gas_pos,
    mass=gas_mass,
    r="infinity",
    xsize=len(x_edges)-1,
    ysize=len(y_edges)-1,
    p=0,  # Viewing angle theta
    roll=0,  # Viewing angle phi
    plot=False,
    logscale=False,
    hsml=2*gas_hsml
)

map_mass2 = qv.get_image()
extent_mass = qv.get_extent()

gas_SFR[gas_SFR<=0] = 1e-10

qv = QuickView(
    gas_pos,
    mass=gas_SFR,
    r="infinity",
    xsize=len(x_edges)-1,
    ysize=len(y_edges)-1,
    p=0,  # Viewing angle theta
    roll=0,  # Viewing angle phi
    plot=False,
    logscale=False,
    hsml=2*gas_hsml
)

map_SFR2 = qv.get_image()
extent_SFR = qv.get_extent()

# Mass map 2
plt.figure()
plt.subplot(111)
plt.imshow(np.log10(map_mass2),extent=extent_mass)
plt.colorbar()
plt.xlim(-12, 12)
plt.ylim(-12, 12)
plt.xlabel("${\\rm Pos}~x~[{\\rm kpc}]$", labelpad=0)
plt.ylabel("${\\rm Pos}~y~[{\\rm kpc}]$", labelpad=-3)
plt.savefig("Map_mass_SPHVIEWER_subkpc.png", dpi=200)

plot_map_SFR2 = np.zeros(np.shape(map_SFR2))
plot_map_SFR2[map_SFR2>1e-8] = map_SFR2[map_SFR2>1e-8]
# SF map 2
plt.figure()
plt.subplot(111)
plt.imshow(np.log10(plot_map_SFR2),extent=extent_SFR, vmax = -.5, vmin=-4.5)
plt.colorbar()
plt.xlim(-12, 12)
plt.ylim(-12, 12)
plt.xlabel("${\\rm Pos}~x~[{\\rm kpc}]$", labelpad=0)
plt.ylabel("${\\rm Pos}~y~[{\\rm kpc}]$", labelpad=-3)
plt.savefig("Map_SFR_SPHVIEWER_subkpc.png", dpi=200)

# Make the KS plot using SPHviewer
plot_map_SFR2 *= 1e2
plot_map_SFR2[plot_map_SFR2<=0] = 1e-6

plt.rcParams.update({    "figure.figsize": (3.15, 3.15),     "figure.subplot.left": 0.18,})
plt.figure()
plt.subplot(111, xscale="log", yscale="log")
plt.scatter(map_mass2.flatten() / 1e6*1e2, plot_map_SFR2.flatten(), s=0.4)
plt.plot(Sigma_g,Sigma_star,'--r',linewidth=.5)
plt.xlim(0.3, 900)
plt.ylim(3e-7, 3)
plt.xlabel("$\\Sigma_g~[{\\rm M_\\odot\\cdot pc^{-2}}]$", labelpad=0)
plt.ylabel("$\\Sigma_{\\rm SFR}~[{\\rm M_\\odot \\cdot yr^{-1} \\cdot kpc^{-2}}]$", labelpad=0)
plt.savefig("KS_law_SPHVIEWER_subkpc.png", dpi=200)


data_mass = 1e-12*np.ones(np.shape(map_mass2))
data_mass[map_mass2>0] = map_mass2[map_mass2>0]
logmap_mass = np.log10(data_mass / 1e6*1e2)
logmap_SFR = np.log10(plot_map_SFR2.flatten())

plt.rcParams.update({    "figure.figsize": (3.15, 3.15),     "figure.subplot.left": 0.15,})
plt.figure()
plt.subplot(111)
plt.hist2d(np.log10(map_mass2.flatten() / 1e6*1e2), np.log10(plot_map_SFR2.flatten()),bins=50,range=[[-1,4],[-5,2]], norm=LogNorm())
plt.plot(np.log10(Sigma_g),np.log10(Sigma_star),'--r',linewidth=2,label="KS law")
plt.plot([-1,4],[-4,1],":b",label="Const. gas depletion time")
plt.xlabel("$\\log \\Sigma_g~[{\\rm M_\\odot\\cdot pc^{-2}}]$", labelpad=0)
plt.ylabel("$\\log \\Sigma_{\\rm SFR}~[{\\rm M_\\odot \\cdot yr^{-1} \\cdot kpc^{-2}}]$", labelpad=0)
plt.legend()
plt.colorbar()
plt.savefig("KS_law_density.png",dpi=300)
plt.close()
