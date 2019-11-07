#!/usr/bin/env python3
import numpy as np
import h5py as h5
import matplotlib.pyplot as plt
import sys

name = str(sys.argv[1])
data = np.loadtxt("SNIa.txt")

g = h5.File("output_%.4d.hdf5" % 1, "r")

unit_time_in_cgs = g["/Units"].attrs["Unit time in cgs (U_t)"]

h = g["/PartType4/SmoothingLengths"][:]
numbstars = len(h)
hmaximalstars = len(h[h == 10.0])
print(hmaximalstars, numbstars)

Myr_in_cgs = 3.154e13

times = data[:, 2] * unit_time_in_cgs / Myr_in_cgs

nu = float(g["/Parameters"].attrs["SNIaDTD:SNIa_efficiency_p_Msun"])
DTD_delay = float(g["/Parameters"].attrs["SNIaDTD:SNIa_delay_time_Gyr"])
DTD_norm = float(g["/Parameters"].attrs["SNIaDTD:normalization_timescale_Gyr"])

StellarMass = (
    np.sum(g["/PartType4/Masses"][:]) - np.sum(g["/PartType4/Masses"][h == 10.0])
) * 1e10
times_expected = np.linspace(0, times[-1], 100)
SNIas = np.zeros(len(times_expected))

for i in range(1, len(SNIas)):
    # We shifted the age of the stars by 97.78 Myr so we need to add this in the prediction calculation
    SNIas[i] = (
        StellarMass
        * nu
        / np.log(DTD_norm / DTD_delay)
        * np.log((times_expected[i] + 97.78) / (times_expected[i - 1] + 97.78))
    )

SNIas_sim = data[:, 9]
SNIas_number = data[:, 12]
SNIas_number_cum = np.cumsum(SNIas_number)
deviation_frac = np.cumsum(SNIas_number) ** 0.5 / np.cumsum(SNIas_number)

plt.plot(times_expected, np.cumsum(SNIas))
plt.plot(times, np.cumsum(SNIas_sim))
plt.fill_between(
    times,
    np.cumsum(SNIas_sim) * (1.0 - deviation_frac),
    np.cumsum(SNIas_sim) * (1.0 + deviation_frac),
    facecolor="#ff7f0e",
    alpha=0.25,
)
plt.xlabel("Time (Myr)")
plt.ylabel("Number of SNIa")
plt.savefig("./SNIa_rate.png")
