###############################################################################
# This file is part of SWIFT.
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

import matplotlib

matplotlib.use("Agg")
from pylab import *
from scipy import stats
import h5py
import numpy as np
import glob
import os.path

# Plot parameters
params = {
    "axes.labelsize": 10,
    "axes.titlesize": 10,
    "font.size": 12,
    "legend.fontsize": 12,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "text.usetex": True,
    "figure.figsize": (9.90, 6.45),
    "figure.subplot.left": 0.1,
    "figure.subplot.right": 0.99,
    "figure.subplot.bottom": 0.1,
    "figure.subplot.top": 0.95,
    "figure.subplot.wspace": 0.2,
    "figure.subplot.hspace": 0.2,
    "lines.markersize": 6,
    "lines.linewidth": 3.0,
    "text.latex.unicode": True,
}
rcParams.update(params)
rc("font", **{"family": "sans-serif", "sans-serif": ["Times"]})

# Number of snapshots and elements
newest_snap_name = max(glob.glob("output_*.hdf5"), key=os.path.getctime)
n_snapshots = int(newest_snap_name.replace("output_", "").replace(".hdf5", "")) + 1

# Read final snapshot
i = n_snapshots-1
print("reading snapshot " + str(i))
# Read the simulation data
sim = h5py.File("output_%04d.hdf5" % i, "r")
star_abundances = sim["/PartType4/ElementMassFractions"][:][:]

PhysicalConstants = sim["/PhysicalConstants/CGS/"]
mp_in_cgs = float(PhysicalConstants.attrs["proton_mass"])

mH_in_cgs = 1.00784*mp_in_cgs
mEu_in_cgs = 151.964*mp_in_cgs
mFe_in_cgs = 55.845*mp_in_cgs
mO_in_cgs = 15.999*mp_in_cgs

# Asplund et al. (2009)
Fe_H_Sun = 7.5
O_H_Sun = 8.69
Eu_H_Sun = 0.52
O_Fe_Sun = O_H_Sun-Fe_H_Sun-np.log10(mFe_in_cgs/mO_in_cgs)
Eu_Fe_Sun = Eu_H_Sun-Fe_H_Sun-np.log10(mFe_in_cgs/mEu_in_cgs)
Fe_H_Sun = Fe_H_Sun-12.0-np.log10(mH_in_cgs/mFe_in_cgs)

Fe_H = np.log10(star_abundances[:,8]/star_abundances[:,0])-Fe_H_Sun
O_Fe = np.log10(star_abundances[:,4]/star_abundances[:,8])-O_Fe_Sun
Eu = star_abundances[:,9]
Eu[Eu==0.] = np.min(Eu[Eu>0]) #Make lower limit
Eu_Fe = np.log10(Eu/star_abundances[:,8])-Eu_Fe_Sun

# Plot the interesting quantities
figure()

# Box stellar abundance --------------------------------
subplot(221)
grid(True)
plot(Fe_H, Eu_Fe,'o',ms=1.5,color='royalblue')
xlabel("[Fe/H]", labelpad=0)
ylabel("[Eu/Fe]", labelpad=0)

# Box stellar abundance --------------------------------
subplot(222)
grid(True)
plot(Fe_H, O_Fe,'o',ms=1.5,color='royalblue')
xlabel("[Fe/H]", labelpad=0)
ylabel("[O/Fe]", labelpad=0)

# -------------------

file = './r_processes.txt'
data = np.loadtxt(file)
time = data[:,2]*9.778131e+02 #Myr
injected_Mass = data[:,8]*1.000347e+10 #Msu
rate = data[:,12]*1.022690e-03 #Num/Myr

# Box event rates --------------------------------------
subplot(223)
grid(True)
plot(time,rate,'-',lw=2,color="crimson")
xlabel("Time [Myr]", labelpad=0)
ylabel("Number of r-processes per time [Myr$^{-1}$]", labelpad=2)
ticklabel_format(style="sci", axis="y", scilimits=(0, 0))

# Box mass --------------------------------
subplot(224)
grid(True)
plot(time,np.cumsum(injected_Mass),'-',lw=2,color="crimson")
xlabel("Time [Myr]", labelpad=0)
ylabel("Cumulative injected mass of Europium [M$_{\odot}$]", labelpad=2)
ticklabel_format(style="sci", axis="y", scilimits=(0, 0))


savefig("Europium_enrichment.png", dpi=200)
