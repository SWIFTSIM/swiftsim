#!/usr/bin/env python3
"""
Read a logger file by using an index file.
Example: ./reader_example.py -t 0.1 ../../examples/SedovBlast_3D/index_*dump
"""
import sys
import numpy as np
import matplotlib.pyplot as plt
from glob import glob
from mpl_toolkits.mplot3d import Axes3D
import argparse
sys.path.append("../.libs/")

import liblogger as logger


def parse_arguments():
    parser = argparse.ArgumentParser(
        description='Read a logfile and plots some basic properties')

    default_files = "../../examples/HydroTests/SedovBlast_3D/index_*dump"
    default_files = glob(default_files)

    parser.add_argument("-t", '--time', dest='time',
                        type=float, default=0.01,
                        help='Simulation time')
    parser.add_argument('files', metavar='filenames', type=str, nargs="*",
                        help='The filenames of the logfiles')
    args = parser.parse_args()
    if len(args.files) == 0:
        args.files = default_files
    return args

def plot3D(pos):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    ax.plot(pos[:, 0], pos[:, 1], pos[:, 2], ".", markersize=0.1)


def plot2D(pos, entropies):
    plt.figure()
    center = np.array([0.5]*3)
    r2 = np.sum((pos - center)**2, axis=1)

    # plot entropy vs distance
    plt.plot(np.sqrt(r2), entropies, '.',
             markersize=0.2)

    plt.xlabel("Radius")
    plt.ylabel("Entropy")


args = parse_arguments()
print("basename: %s" % args.files)
print("time: %g" % args.time)

# read the logger
positions = np.empty((0, 3))
entropies = np.empty(0)
gas_type = 0
for f in args.files:
    if f.endswith(".dump"):
        filename = f[:-5]
    else:
        raise Exception("It seems that you are not providing a logfile (.dump)")
    with logger.Reader(filename, verbose=0) as reader:


        # Ensure that the fields are present
        fields = ["Coordinates", "Entropies", "ParticleIDs"]
        missing = set(fields).difference(set(reader.get_list_fields(part_type=gas_type)))
        if missing:
            raise Exception("Fields %s not found in the logfile." % missing)

        if ("Coordinates" not in fields or
            "Entropies" not in fields):
            raise Exception("Field not found in the logfile")

        # Read all the particles
        t = reader.get_time_limits()
        out = reader.get_particle_data(
            fields=fields, time=args.time, part_type=gas_type)

        # Get the particle ids
        gas_ids = out[2]
        gas_ids = gas_ids[:len(gas_ids)//2]
        ids = [None] * 6
        ids[gas_type] = gas_ids

        # Read from the ids
        # As we are filtering by particle ids, the field "ParticleIDs" is required
        # in order to verify the particle obtained.
        out = reader.get_particle_data(
            fields=fields, time=args.time, filter_by_ids=ids)

        # Print the missing ids
        gas_ids, ids_found = set(gas_ids), set(out[2])
        print("The following ids were not found: ", gas_ids.difference(ids_found))
        print("The following ids are wrongly missing: ", ids_found.difference(gas_ids))

        # add the data to the list
        positions = np.append(positions, out[0], axis=0)
        entropies = np.append(entropies, out[1], axis=0)

print("Min/Max of the position:", positions.min(), positions.max())
print("Min/Max of the entropy:", entropies.min(), entropies.max())
plot3D(positions)
plot2D(positions, entropies)
plt.show()
