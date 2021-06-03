#!/usr/bin/env python3

# -----------------------------------------------------------------------
# Collection of sanity checks for the 'GEAR' RT scheme in swift for the
# Any run done with the 'GEAR' scheme must pass these tests. Swift must
# be compiled with the '--enable-debugging-checks' flag.
#
# Usage:
#   ./rt_uniform_box_checks-GEAR.py
# or
#   ./rt_uniform_box_checks-GEAR.py snapshot_basename
#
# where 'snapshot_basename' is the base name of the snapshots that you
# want to work with. (Same argument as Snapshots:basename in the .yml
# file). If 'snapshot_basename' is not given, this script assumes it to
# be 'output'.
#
# The script will then read all the snapshot_basename_XXXX.hdf5 files
# in the directory where you're running the script from
# -----------------------------------------------------------------------


import numpy as np
from sys import argv
from swift_rt_GEAR_io import get_snap_data


# some behaviour options
print_diffs = False  # print differences you find
print_additional_information = False
break_on_diff = False  # quit when you find a difference

# tolerance for a float to be equal
float_comparison_tolerance = 1e-4


if len(argv) > 1:
    file_prefix = argv[1]
else:
    file_prefix = "output"


def check_injection(snapdata, rundata):
    """
    Do checks related to energy injections.

    snapdata: list of swift_rt_GEAR_io.RTSnapData objects

    rundata: swift_rt_GEAR_io.Rundata object
    """

    print("checking injection")

    # ----------------------------------------------------------------
    # Check 1: Make sure the right amount of energy has been injected
    # into the gas
    # ----------------------------------------------------------------

    initial_energies = np.sum(snapdata[0].gas.PhotonEnergies, axis=1)
    initial_time = snapdata[0].time

    emission_rates = rundata.const_emission_rates
    ngroups = rundata.ngroups

    # Check 1a) : sum initial energy + sum injected energy = sum current energy
    # --------------------------------------------------------------------------
    for snap in snapdata[1:]:
        dt = snap.time - initial_time
        photon_energies = np.sum(snap.gas.PhotonEnergies, axis=1)
        injected_energies = np.sum(snap.stars.InjectedPhotonEnergy, axis=0)

        for g in range(ngroups):
            energy_expected = initial_energies[g] + photon_energies[g]
            diff = abs(1.0 - photon_energies[g] / energy_expected)
            if diff > float_comparison_tolerance:
                print(
                    "--- Injection Energy Budget is wrong; snapshot",
                    snap.snapnr,
                    "group",
                    g + 1,
                    "diff",
                    diff,
                    photon_energies[g],
                    energy_expected,
                )
                if break_on_diff:
                    quit()

    # Check 1b) : sum injected energy >= (t_now - t_start * injection_rate)
    # --------------------------------------------------------------------------
    # we may have injected too much energy, because stars inject all the
    # radiation of their entire time step instantaneously

    # TODO: this assumes a constant number of stars

    if rundata.use_const_emission_rate and not rundata.hydro_controlled_injection:
        if snapdata[0].snapnr != 0:
            print("You need to read in snapshot 0 to do this particular test")
        else:

            for snap in snapdata[1:]:
                dt = snap.time - initial_time
                injected_energies = np.sum(snap.stars.InjectedPhotonEnergy, axis=0)
                energies_expected = snap.nstars * emission_rates * dt
                energies_expected = energies_expected.to(injected_energies.units)

                if print_additional_information:
                    print(
                        "Snapshot",
                        snap.snapnr,
                        "Overshoot in injection due to discretization:",
                        (injected_energies / energies_expected).value,
                        "[This is expected behaviour]",
                    )

                for g in range(ngroups):
                    if energies_expected[g] > injected_energies[g]:
                        diff = abs(1.0 - injected_energies[g] / energies_expected[g])
                        if diff > float_comparison_tolerance:
                            print("--- Injected Energy Prediction is wrong;")
                            print(
                                "--- Snapshot",
                                snap.snapnr,
                                "group",
                                g + 1,
                                "injected:",
                                injected_energies[g],
                                "expected:",
                                energies_expected[g],
                                "ratio",
                                injected_energies[g] / energies_expected[g],
                                "diff",
                                diff,
                            )
                            if break_on_diff:
                                quit()


def main():
    """
    Main function to run.
    """

    snapdata, rundata = get_snap_data(prefix=file_prefix)

    check_injection(snapdata, rundata)

    return


if __name__ == "__main__":
    main()
