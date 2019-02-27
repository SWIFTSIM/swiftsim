"""
Plots the solution for the ContactDiscontinuity_1D test.
"""

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from swiftsimio import load

matplotlib.use("Agg")


def get_data_dump(metadata):
    """
    Gets a big data dump from the SWIFT metadata.
    """

    output = (
        "SWIFT\n"
        + metadata.code_info
        + "\n\n"
        + "Compiler\n"
        + metadata.compiler_info
        + "\n\n"
        + "Hydrodynamics\n"
        + metadata.hydro_info
        + "\n\n"
        + "Viscosity\n"
        + metadata.viscosity_info
    )

    return output


def get_data_list(start: int, stop: int, handle: str):
    """
    Gets a list of swiftsimio objects that contains all of the data.
    """

    data = [load("{}_{:04d}.hdf5".format(handle, x)) for x in range(start, stop + 1)]

    return data


def setup_axes(size=[8, 4], dpi=300):
    """
    Sets up the axes with the correct labels, etc.
    """
    fig, ax = plt.subplots(ncols=2, figsize=size, dpi=dpi)

    ax.flatten()

    ax[0].axis("off")
    ax[1].set_xlabel("Time [s]")
    ax[1].set_ylabel("Relative energy difference $u/\\left<u\\right>$")

    return fig, ax


def extract_plottables(data_list):
    """
    Extracts the plottables. Returns:
    mean, stdev, max, min
    """

    data = [
        np.diff(x.gas.internal_energy.value) / np.mean(x.gas.internal_energy.value)
        for x in data_list
    ]

    means = np.array([np.mean(x) for x in data])
    stdevs = np.array([np.std(x) for x in data])
    maxs = np.array([np.max(x) for x in data])
    mins = np.array([np.min(x) for x in data])

    return means, stdevs, maxs, mins


def make_plot(start: int, stop: int, handle: str):
    """
    Makes the plot and returns the figure and axes objects.
    """
    fig, ax = setup_axes()
    data_list = get_data_list(start, stop, handle)
    data_dump = get_data_dump(data_list[0].metadata)
    t = [x.metadata.t for x in data_list]
    means, stdevs, maxs, mins = extract_plottables(data_list)

    ax[0].text(
        0.5,
        0.5,
        data_dump,
        ha="center",
        va="center",
        fontsize=8,
        transform=ax[0].transAxes,
    )

    ax[1].fill_between(
        t, means - stdevs, means + stdevs, color="C0", alpha=0.5, edgecolor="none"
    )
    ax[1].plot(t, means, label="Mean", c="C0")
    ax[1].plot(t, maxs, label="Max", linestyle="dashed", c="C1")
    ax[1].plot(t, mins, label="Min", linestyle="dashed", c="C2")

    ax[1].legend(loc=1, markerfirst=False)

    ax[1].set_xlim(t[0], t[-1])

    fig.tight_layout()

    return fig, ax


if __name__ == "__main__":
    import argparse as ap

    parser = ap.ArgumentParser(
        description="Makes a plot of the data from the ContactDiscontinuity_1D test."
    )

    parser.add_argument(
        "-i", "--initial", help="Initial snapshot. Default: 0", type=int, default=0
    )

    parser.add_argument(
        "-f", "--final", help="Final snapshot. Default: 50", type=int, default=50
    )

    parser.add_argument(
        "-s",
        "--snapshot",
        help="First part of the snapshot filename. Default: contact",
        type=str,
        default="contact",
    )

    parser.add_argument(
        "-o",
        "--output",
        help="Output filename. Default: contact.png",
        type=str,
        default="contact.png",
    )

    args = vars(parser.parse_args())

    fig, ax = make_plot(args["initial"], args["final"], args["snapshot"])

    fig.savefig(args["output"])
