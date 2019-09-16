#!/usr/bin/env python3
import numpy as np
import h5py as h5
import matplotlib.pyplot as plt
import sys

name = str(sys.argv[1])
data = np.loadtxt("SNIa.txt")

g = h5.File("output_%.4d.hdf5" % 1, "r")

h = g["/PartType4/SmoothingLengths"][:]
numbstars = len(h)
hmaximalstars = len(h[h == 10.0])
print(hmaximalstars, numbstars)

times = data[:, 2] * 9.778131e+02

nu = float(g["/Parameters"].attrs["SNIaDTD:SNIa_efficiency_p_Msun"])
DTD_delay = float(g["/Parameters"].attrs["SNIaDTD:SNIa_delay_time_Gyr"])
DTD_norm = float(g["/Parameters"].attrs["SNIaDTD:normalization_timescale_Gyr"])

StellarMass = (np.sum(g["/PartType4/Masses"][:])-np.sum(g["/PartType4/Masses"][h == 10.0]))*1e10
times_expected = np.linspace(0,times[-1],100)
SNIas = np.zeros(len(times_expected))

for i in range(1,len(SNIas)):
    SNIas[i] = StellarMass * nu / np.log(DTD_norm/DTD_delay) * np.log((times_expected[i]+97.78)/(times_expected[i-1]+97.78))

SNIas_sim = data[:,9]
SNIas_number = data[:,12]
SNIas_number_cum = np.cumsum(SNIas_number)
deviation_frac = np.cumsum(SNIas_number)**.5 /np.cumsum(SNIas_number) 

plt.plot(times_expected, np.cumsum(SNIas))
plt.plot(times, np.cumsum(SNIas_sim))
plt.fill_between(
        times,
        np.cumsum(SNIas_sim)*(1.-deviation_frac),
        np.cumsum(SNIas_sim)*(1.+deviation_frac),
        facecolor="#ff7f0e",
        alpha=0.25)
plt.xlabel("Time (Myr)")
plt.ylabel("Number of SNIa")
plt.savefig("./SNIa_rate.png")
