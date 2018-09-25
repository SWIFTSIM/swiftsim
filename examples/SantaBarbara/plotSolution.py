"""
Plots the "solution" (i.e. some profiles) for the Santa Barbara cluster.

Invoke as follows:

python3 plotSolution.py <snapshot number> <catalogue directory> <number of bins (optional)>
"""

import matplotlib.pyplot as plt
import numpy as np

import h5py

from collections import namedtuple
from typing import Tuple

# Simulation data
SimulationData = namedtuple("SimulationData", ["gas", "dark_matter", "metadata"])
ParticleData = namedtuple(
    "ParticleData", ["coordinates", "radii", "masses", "densities", "energies"]
)
MetaData = namedtuple("MetaData", ["header", "code", "hydroscheme"])


def get_energies(handle: h5py.File):
    """
    Gets the energies with the correct units.
    """
    u = handle["PartType0/InternalEnergy"][:]
    unit_length_in_cgs = handle["/Units"].attrs["Unit length in cgs (U_L)"]
    unit_mass_in_cgs = handle["/Units"].attrs["Unit mass in cgs (U_M)"]
    unit_time_in_cgs = handle["/Units"].attrs["Unit time in cgs (U_t)"]
    gas_gamma = handle["/HydroScheme"].attrs["Adiabatic index"][0]
    a = handle["/Cosmology"].attrs["Scale-factor"][0]

    unit_length_in_si = 0.01 * unit_length_in_cgs
    unit_mass_in_si = 0.001 * unit_mass_in_cgs
    unit_time_in_si = unit_time_in_cgs

    u *= unit_length_in_si ** 2 / unit_time_in_si ** 2
    u /= a ** (3 * (gas_gamma - 1.))

    return u


def load_data(filename: str, center: np.array) -> SimulationData:
    """
    Loads the relevant data for making the profiles, as well as some metadata
    for the plot.

    Center is the center of the SB cluster and is used to calculate the radial
    distances to the particles.
    """

    with h5py.File(filename, "r") as file:
        gas_handle = file["PartType0"]
        dm_handle = file["PartType1"]

        gas_data = ParticleData(
            coordinates=gas_handle["Coordinates"][...],
            radii=get_radial_distances(gas_handle["Coordinates"][...], center),
            masses=gas_handle["Masses"][...],
            energies=get_energies(file),
            densities=gas_handle["Density"][...],
        )

        dm_data = ParticleData(
            coordinates=dm_handle["Coordinates"][...],
            radii=get_radial_distances(dm_handle["Coordinates"][...], center),
            masses=gas_handle["Masses"][...],
            energies=None,
            densities=None,
        )

        metadata = MetaData(
            header=dict(file["Header"].attrs),
            code=dict(file["Code"].attrs),
            hydroscheme=dict(file["HydroScheme"].attrs),
        )

        simulation_data = SimulationData(
            gas=gas_data, dark_matter=dm_data, metadata=metadata
        )

    return simulation_data


def get_halo_center(catalogue_filename: str) -> np.array:
    """
    Gets the halo center of the largest halo (i.e. the SB cluster).

    You will want the .properties file, probably

    halo/santabarbara.properties

    that is given by VELOCIraptor.
    """

    with h5py.File(catalogue_filename, "r") as file:
        x = file["Xc"][0]
        y = file["Yc"][0]
        z = file["Zc"][0]

    return np.array([x, y, z])


def get_radial_distances(coordinates: np.ndarray, center: np.array) -> np.array:
    """
    Gets the radial distances for all particles.
    """
    dx = coordinates - center

    return np.sqrt(np.sum(dx * dx, axis=1))


def get_radial_density_profile(radii, masses, bins: int) -> Tuple[np.ndarray]:
    """
    Gets the radial gas density profile, after generating similar bins to those
    used in similar works.
    """

    bins = np.logspace(-2, 1, bins)

    histogram, bin_edges = np.histogram(a=radii, weights=masses, bins=bins)

    volumes = np.array(
        [
            (4. * np.pi / 3.) * (r_outer ** 3 - r_inner ** 3)
            for r_outer, r_inner in zip(bin_edges[1:], bin_edges[:-1])
        ]
    )

    return histogram / volumes, bin_edges  # densities


def mu(T, H_frac, T_trans):
    """
    Get the molecular weight as a function of temperature.
    """
    if T > T_trans:
        return 4. / (8. - 5. * (1. - H_frac))
    else:
        return 4. / (1. + 3. * H_frac)


def T(u, metadata: MetaData):
    """
    Temperature of primordial gas.
    """

    gas_gamma = metadata.hydroscheme["Adiabatic index"][0]
    H_frac = metadata.hydroscheme["Hydrogen mass fraction"][0]
    T_trans = metadata.hydroscheme["Hydrogen ionization transition temperature"][0]

    k_in_J_K = 1.38064852e-23
    mH_in_kg = 1.6737236e-27

    T_over_mu = (gas_gamma - 1.) * u * mH_in_kg / k_in_J_K
    ret = np.ones(np.size(u)) * T_trans

    # Enough energy to be ionized?
    mask_ionized = T_over_mu > (T_trans + 1) / mu(T_trans + 1, H_frac, T_trans)
    if np.sum(mask_ionized) > 0:
        ret[mask_ionized] = T_over_mu[mask_ionized] * mu(T_trans * 10, H_frac, T_trans)

    # Neutral gas?
    mask_neutral = T_over_mu < (T_trans - 1) / mu((T_trans - 1), H_frac, T_trans)
    if np.sum(mask_neutral) > 0:
        ret[mask_neutral] = T_over_mu[mask_neutral] * mu(0, H_frac, T_trans)

    return ret


def get_radial_temperature_profile(data: SimulationData, bins: int) -> np.ndarray:
    """
    Gets the radial gas density profile, after generating similar bins to those
    used in similar works.
    """

    temperatures = T(data.gas.energies, data.metadata)
    radii = data.gas.radii

    bins = np.logspace(-2, 1, bins)

    histogram, _ = np.histogram(a=radii, weights=temperatures, bins=bins)

    counts, _ = np.histogram(a=radii, weights=np.ones_like(radii), bins=bins)

    return histogram / counts  # need to get mean value in bin


def get_radial_entropy_profile(data: SimulationData, bins: int) -> np.ndarray:
    """
    Gets the radial gas density profile, after generating similar bins to those
    used in similar works.
    """

    gas_gamma = data.metadata.hydroscheme["Adiabatic index"][0]
    gamma_minus_one = gas_gamma - 1.0

    entropies = (
        data.gas.energies * (gamma_minus_one) / data.gas.densities ** gamma_minus_one
    )
    print("Warning: Current entropy profile assumes all gas is ionised")
    radii = data.gas.radii

    bins = np.logspace(-2, 1, bins)

    histogram, _ = np.histogram(a=radii, weights=entropies, bins=bins)

    counts, _ = np.histogram(a=radii, weights=np.ones_like(radii), bins=bins)

    return histogram / counts  # need to get mean value in bin


def create_plot(data: SimulationData, bins: int):
    """
    Creates the figure and axes objects and plots the data on them.
    """

    fig, axes = plt.subplots(2, 2, figsize=(8, 8), sharex=True)
    axes = axes.flatten()

    gas_density, bin_edges = get_radial_density_profile(
        data.gas.radii, data.gas.masses, bins=bins
    )
    dm_density, _ = get_radial_density_profile(
        data.dark_matter.radii, data.dark_matter.masses, bins=bins
    )
    temperature = get_radial_temperature_profile(data, bins=bins)
    entropy = get_radial_entropy_profile(data, bins=bins)

    bin_centers = [0.5 * (x + y) for x, y in zip(bin_edges[:-1], bin_edges[1:])]

    axes[0].loglog()
    axes[0].scatter(bin_centers, gas_density)
    axes[0].set_ylabel(r"$\rho_{\rm gas} (R)$ [$10^{10}$ M$_\odot$ Mpc$^{-3}$]")
    axes[0].set_xlabel(r"R [Mpc]")

    axes[1].scatter(bin_centers, entropy / 10000)
    axes[1].set_ylabel(r"Entropy $A(R)$ [Arbirtary]")
    axes[1].set_xlabel(r"R [Mpc]")

    axes[2].loglog()
    axes[2].scatter(bin_centers, temperature)
    axes[2].set_ylabel(r"$T_{\rm gas} (R)$ [K]")
    axes[2].set_xlabel(r"R [Mpc]")

    axes[3].loglog()
    axes[3].scatter(bin_centers, dm_density)
    axes[3].set_ylabel(r"$\rho_{\rm DM} (R)$ [$10^{10}$ M$_\odot$ Mpc$^{-3}$]")
    axes[3].set_xlabel(r"R [Mpc]")

    fig.tight_layout()

    return fig, axes


if __name__ == "__main__":
    import sys

    filename = "santabarbara_{:04d}.hdf5".format(int(sys.argv[1]))
    catalogue_filename = f"{sys.argv[2]}/santabarbara.properties"

    try:
        bins = int(sys.argv[3])
    except:
        bins = 25

    halo_center = get_halo_center(catalogue_filename)
    simulation_data = load_data(filename, halo_center)

    fig, _ = create_plot(data=simulation_data, bins=bins)

    fig.savefig("santabarbara.png", dpi=300)
