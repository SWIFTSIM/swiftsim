#!/usr/bin/env python3
################################################################################
# This file is part of SWIFT.
# Copyright (c) 2019 Folkert Nobels    (nobels@strw.leidenuniv.nl)
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
import numpy as np
import h5py as h5

import matplotlib

matplotlib.use("Agg")
font = {"size": 6}
matplotlib.rc("font", **font)
import matplotlib.pyplot as plt
from sphviewer.tools import QuickView, cmaps
from matplotlib.colors import LogNorm
from matplotlib import gridspec
import sphviewer as sph
import sys

# snapshot
snap = int(sys.argv[1])
# size in kpc of radius of interest
bs = float(sys.argv[2])

filename = "output_%.4d.hdf5" % snap
fid = "output_%.4d.hdf5" % 0

f = h5.File(filename, "r")
g = h5.File(fid, "r")

# Physical constants
k_in_cgs = 1.38064852e-16
mH_in_cgs = 1.6737236e-24
year_in_cgs = 3600.0 * 24 * 365.0
Msun_in_cgs = 1.98848e33
G_in_cgs = 6.67259e-8
pc_in_cgs = 3.08567758e18
Msun_p_pc2 = Msun_in_cgs / pc_in_cgs ** 2

# Read units
unit_length_in_cgs = f["/Units"].attrs["Unit length in cgs (U_L)"]
unit_mass_in_cgs = f["/Units"].attrs["Unit mass in cgs (U_M)"]
unit_time_in_cgs = f["/Units"].attrs["Unit time in cgs (U_t)"]

boxsize = f["/Header"].attrs["BoxSize"]
centre = boxsize / 2.0

# Read gas properties
gas_pos = f["/PartType0/Coordinates"][:, :]
gas_mass = f["/PartType0/Masses"][:]
gas_rho = f["/PartType0/Densities"][:]
gas_T = f["/PartType0/Temperatures"][:]
gas_hsml = f["/PartType0/SmoothingLengths"][:]
gas_v = f["/PartType0/Velocities"][:, :]
gas_XH = f["/PartType0/ElementMassFractions"][:, 0]
gas_ID = f["/PartType0/ParticleIDs"][:]

# Read the 0 particles
gas_pos0 = g["/PartType0/Coordinates"][:, :]
gas_mass0 = g["/PartType0/Masses"][:]
gas_rho0 = g["/PartType0/Densities"][:]
gas_T0 = g["/PartType0/Temperatures"][:]
gas_hsml0 = g["/PartType0/SmoothingLengths"][:]
gas_v0 = g["/PartType0/Velocities"][:, :]

# Turn the mass into better units
gas_mass *= unit_mass_in_cgs / Msun_in_cgs
gas_mass0 *= unit_mass_in_cgs / Msun_in_cgs

# Centre the box
gas_pos[:, 0] -= centre[0]
gas_pos[:, 1] -= centre[1]
gas_pos[:, 2] -= centre[2]
gas_pos0[:, 0] -= centre[0]
gas_pos0[:, 1] -= centre[1]
gas_pos0[:, 2] -= centre[2]

gas_rho = gas_rho * unit_mass_in_cgs / unit_length_in_cgs ** 3
gas_rho0 = gas_rho0 * unit_mass_in_cgs / unit_length_in_cgs ** 3
gas_nH = gas_rho  # * unit_mass_in_cgs / unit_length_in_cgs ** 3
gas_nH /= mH_in_cgs
gas_nH *= gas_XH

intervals = 501
Particles = sph.Particles(gas_pos, gas_mass, 2 * gas_hsml)
Camera = sph.Camera(
    r="infinity",
    t=0,
    p=0,
    roll=0,
    xsize=intervals,
    ysize=intervals,
    x=0,
    y=0,
    z=0,
    extent=[-bs, bs, -bs, bs],
)
Scene = sph.Scene(Particles, Camera)
Render = sph.Render(Scene)

extent = Render.get_extent()
lpixel = (extent[1] - extent[0]) / intervals

density_image = Render.get_image()

# fiducial positions
x = np.linspace(-bs, bs, intervals)

# Temperature phase on

intervals = 501
Particles = sph.Particles(gas_pos, gas_T * gas_mass, 2 * gas_hsml)
Camera = sph.Camera(
    r="infinity",
    t=0,
    p=0,
    roll=0,
    xsize=intervals,
    ysize=intervals,
    x=0,
    y=0,
    z=0,
    extent=[-bs, bs, -bs, bs],
)
Scene = sph.Scene(Particles, Camera)
Render = sph.Render(Scene)

extent = Render.get_extent()
lpixel = (extent[1] - extent[0]) / intervals

T_image = Render.get_image() / density_image

##################

Particles = sph.Particles(gas_pos, gas_mass, 2 * gas_hsml)
Camera = sph.Camera(
    r="infinity",
    t=0,
    p=90,
    roll=90,
    xsize=intervals,
    ysize=intervals,
    x=0,
    y=0,
    z=0,
    extent=[-bs, bs, -bs, bs],
)
Scene = sph.Scene(Particles, Camera)
Render = sph.Render(Scene)

extent = Render.get_extent()

lpixel = (extent[1] - extent[0]) / intervals

density_image2 = Render.get_image()

Particles = sph.Particles(gas_pos, gas_T * gas_mass, 2 * gas_hsml)
Camera = sph.Camera(
    r="infinity",
    t=0,
    p=90,
    roll=90,
    xsize=intervals,
    ysize=intervals,
    x=0,
    y=0,
    z=0,
    extent=[-bs, bs, -bs, bs],
)
Scene = sph.Scene(Particles, Camera)
Render = sph.Render(Scene)

T_image2 = Render.get_image() / density_image2

# fiducial positions
x = np.linspace(-bs, bs, intervals)

padding = -14
fig = plt.figure(figsize=(4.0, 4.0))

ax1 = fig.add_axes([0.05, 0.1, 0.4, 0.4])
i = ax1.imshow(
    np.log10(density_image2) - np.log10(lpixel ** 2) - 6.0,
    extent=[
        extent[0] + Scene.Camera.get_params()["x"],
        extent[1] + Scene.Camera.get_params()["x"],
        extent[2] + Scene.Camera.get_params()["y"],
        extent[3] + Scene.Camera.get_params()["y"],
    ],
    vmin=-3.0,
    vmax=3.0,
)
ax1.set_xticklabels([])
ax1.set_yticklabels([])
plt.tick_params(axis="y", which="both", left=False, right=False)
plt.tick_params(axis="x", which="both", top=False, bottom=False)
colorbar_ax = fig.add_axes([0.45, 0.1, 0.05, 0.4])
cbar = fig.colorbar(i, cax=colorbar_ax)
# cbar.ax.set_yticklabels([-2,-1,0.,1.,2.,3.,4.,5.,6.])
cbar.ax.yaxis.set_tick_params(pad=padding)
cbytick_obj = plt.getp(cbar.ax.axes, "yticklabels")
plt.setp(cbytick_obj, color="w")
yticks = cbar.ax.yaxis.get_major_ticks()
yticks[-1].set_visible(False)
yticks[0].set_visible(False)

ax2 = fig.add_axes([0.05, 0.5, 0.4, 0.4])
i = ax2.imshow(
    np.log10(density_image) - np.log10(lpixel ** 2) - 6.0,
    extent=[
        extent[0] + Scene.Camera.get_params()["x"],
        extent[1] + Scene.Camera.get_params()["x"],
        extent[2] + Scene.Camera.get_params()["y"],
        extent[3] + Scene.Camera.get_params()["y"],
    ],
    vmin=-3,
    vmax=3.0,
)
ax2.set_xticklabels([])
ax2.set_yticklabels([])
plt.hlines(-15, -15, -5, color="w", linewidth=1.5)
plt.text(-10, -18, "10 kpc", color="w", horizontalalignment="center")
plt.text(
    0,
    22,
    "Density $[{\\rm M_\\odot\\cdot pc^{-2}}]$",
    color="w",
    horizontalalignment="center",
)
plt.tick_params(axis="y", which="both", left=False, right=False)
plt.tick_params(axis="x", which="both", top=False, bottom=False)
colorbar_ax = fig.add_axes([0.45, 0.5, 0.05, 0.4])
cbar = fig.colorbar(i, cax=colorbar_ax)
# cbar.ax.set_yticklabels([-2,-1,0.,1.,2.])
cbar.ax.yaxis.set_tick_params(pad=padding)
cbytick_obj = plt.getp(cbar.ax.axes, "yticklabels")
plt.setp(cbytick_obj, color="w")
yticks = cbar.ax.yaxis.get_major_ticks()
yticks[-1].set_visible(False)
yticks[0].set_visible(False)

ax3 = fig.add_axes([0.5, 0.1, 0.4, 0.4])
i = ax3.imshow(
    np.log10(T_image2),
    extent=[
        extent[0] + Scene.Camera.get_params()["x"],
        extent[1] + Scene.Camera.get_params()["x"],
        extent[2] + Scene.Camera.get_params()["y"],
        extent[3] + Scene.Camera.get_params()["y"],
    ],
    vmin=4.0,
    vmax=7.5,
    cmap="plasma",
)
ax3.set_xticklabels([])
ax3.set_yticklabels([])
plt.tick_params(axis="y", which="both", left=False, right=False)
plt.tick_params(axis="x", which="both", top=False, bottom=False)
colorbar_ax = fig.add_axes([0.90, 0.1, 0.05, 0.4])
cbar = fig.colorbar(i, cax=colorbar_ax)
# cbar.ax.set_yticklabels([4.5,5.0,5.5,6.0,7.0])
cbar.ax.yaxis.set_tick_params(pad=padding)
cbytick_obj = plt.getp(cbar.ax.axes, "yticklabels")
plt.setp(cbytick_obj, color="w")
yticks = cbar.ax.yaxis.get_major_ticks()
yticks[-1].set_visible(False)
yticks[0].set_visible(False)

ax4 = fig.add_axes([0.5, 0.5, 0.4, 0.4])
i = ax4.imshow(
    np.log10(T_image),
    extent=[
        extent[0] + Scene.Camera.get_params()["x"],
        extent[1] + Scene.Camera.get_params()["x"],
        extent[2] + Scene.Camera.get_params()["y"],
        extent[3] + Scene.Camera.get_params()["y"],
    ],
    vmin=4.0,
    vmax=7.5,
    cmap="plasma",
)
ax4.set_xticklabels([])
ax4.set_yticklabels([])
plt.hlines(-15, -15, -5, color="w", linewidth=1.5)
plt.text(-10, -18, "10 kpc", color="w", horizontalalignment="center")
plt.text(0, 22, "Temperature [K]", color="w", horizontalalignment="center")
plt.tick_params(axis="y", which="both", left=False, right=False)
plt.tick_params(axis="x", which="both", top=False, bottom=False)
colorbar_ax = fig.add_axes([0.90, 0.5, 0.05, 0.4])
cbar = fig.colorbar(i, cax=colorbar_ax)
# cbar.ax.set_yticklabels([4.5,5.0,5.5,6.0,7.0])
cbar.ax.yaxis.set_tick_params(pad=padding)
cbytick_obj = plt.getp(cbar.ax.axes, "yticklabels")
plt.setp(cbytick_obj, color="w")
yticks = cbar.ax.yaxis.get_major_ticks()
yticks[-1].set_visible(False)
yticks[0].set_visible(False)
# cbar.ax.set_ytick_params(color="w")
plt.savefig("./frames/combined_frame_%.4d.png" % snap, dpi=300)
