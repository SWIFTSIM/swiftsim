#!/usr/bin/env python3

# -----------------------------------------------------------------------
# Collection of checks for the 'GEAR' RT scheme in swift for the
# uniform box test where particles don't move and every time step an
# output file is generated. Swift must be compiled with the
# '--enable-debugging-checks' flag.
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
skip_snap_zero = True  # skip snap_0000.hdf5
skip_last_snap = False  # skip snap_0max.hdf5
skip_coords = True  # skip coordinates check
skip_sml = True  # skip smoothing lengths check
print_diffs = False  # print differences you find
break_on_diff = False  # quit when you find a difference

# tolerance for a float to be equal
float_comparison_tolerance = 1e-5


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

    if not rundata.use_const_emission_rate:
        print("Sim wasn't run with const emission rates. Skipping those tests.")
        return

    # ----------------------------------------------------------------
    # Check 1: Make sure the right amount of energy has been injected
    # into the gas
    # ----------------------------------------------------------------

    initial_energies = [np.sum(s) for s in snapdata[0].gas.PhotonEnergies]
    initial_time = snapdata[0].time

    emission_rates = rundata.const_emission_rates
    ngroups = rundata.ngroups

    for snap in snapdata:
        dt = snap.time - initial_time
        photon_energies = [np.sum(s) for s in snap.gas.PhotonEnergies]

        for g in range(ngroups):
            injected = snap.nstars * emission_rates[g] * dt
            energies_expected = initial_energies[g] + injected
            diff = abs(1.0 - energies_expected / photon_energies[g])
            if diff > float_comparison_tolerance:
                print(
                    "--- Snapshot",
                    snap.snapnr,
                    "group",
                    g + 1,
                    "diff",
                    diff,
                    photon_energies[g],
                    energies_expected,
                )
                if break_on_diff:
                    quit()


def main():
    """
    Main function to run.
    """

    snapdata, rundata = get_snap_data(
        prefix=file_prefix, skip_snap_zero=skip_snap_zero, skip_last_snap=skip_last_snap
    )

    check_injection(snapdata, rundata)

    return


if __name__ == "__main__":
    main()
