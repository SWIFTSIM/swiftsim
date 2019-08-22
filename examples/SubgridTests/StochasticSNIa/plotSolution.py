#!/usr/bin/env python3
import numpy as np
import h5py as h5
import matplotlib.pyplot as plt
import sys

name = str(sys.argv[1])
data = np.loadtxt("SNIa.txt")

g = h5.File("output_%.4d.hdf5" % 1, "r")

h = g["/PartType4/SmoothingLength"][:]
numbstars = len(h)
hmaximalstars = len(h[h == 10.0])

times = data[:, 1] * 0.9785

nu = float(g["/Parameters"].attrs["SNIaDTD:SNIa_efficiency_p_Msun"])
DTD_delay = float(g["/Parameters"].attrs["SNIaDTD:SNIa_delay_time_Gyr"])
DTD_norm = float(g["/Parameters"].attrs["SNIaDTD:normalization_timescale_Gyr"])

N = 1000
N_bins = 1000
t_sim = 0.09785
t_start = 0.09785
time_array = np.linspace(0, t_sim, N) + t_start


def subplotvoid(time_array, t_start, times, numbs, hmaxs, backup=False, backup8=False):
    plt.hist(
        times,
        bins=N_bins,
        cumulative=True,
        histtype="step",
        label="Isolated galaxy simulation",
    )

    # lets plot calculate the ideal curve
    idealcurve = (
        nu / (1.736e-2 * np.log(DTD_norm / DTD_delay)) * np.log((time_array) / t_start)
    ) * (numbs - hmaxs)

    plt.plot(
        time_array - t_start,
        idealcurve,
        label="Ideal prediction power law beta 1",
        color="#ff7f0e",
    )

    plt.fill_between(
        time_array - t_start,
        idealcurve + np.sqrt(idealcurve),
        idealcurve - np.sqrt(idealcurve),
        facecolor="#ff7f0e",
        alpha=0.25,
    )
    plt.fill_between(
        time_array - t_start,
        idealcurve + 2 * np.sqrt(idealcurve),
        idealcurve - 2 * np.sqrt(idealcurve),
        facecolor="#ff7f0e",
        alpha=0.25,
    )
    plt.fill_between(
        time_array - t_start,
        idealcurve + 3 * np.sqrt(idealcurve),
        idealcurve - 3 * np.sqrt(idealcurve),
        facecolor="#ff7f0e",
        alpha=0.25,
    )

    plt.xlabel("Time (Gyr)")
    plt.ylabel("Cumulative number star particles that go SNIa")
    plt.legend()


plt.figure(figsize=(6, 5))
plt.subplot(111)
plt.title("nu = %1.6f, t_delay = %1.6f, t_norm = %1.4f" % (nu, DTD_delay, DTD_norm))
subplotvoid(time_array, t_start, times, numbstars, hmaximalstars, backup=True)
plt.savefig("./" + name + ".png")
