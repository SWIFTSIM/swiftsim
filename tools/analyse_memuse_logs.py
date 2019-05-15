#!/usr/bin/env python
"""
Usage:
    process_memuse.py [options] memuse_report1.dat [memuse_report2.dat] ...

Parse the output of a run of SWIFT to convert the memuse output dumps into a
timeseries of memory use. Also outputs use in memory per labelled type.

This file is part of SWIFT.
Copyright (c) 2019 Peter W. Draper (p.w.draper@durham.ac.uk)

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

from collections import OrderedDict
import argparse
import numpy as np
import sys

#  Command-line arguments.
parser = argparse.ArgumentParser(description="Analyse memory usage reports")

parser.add_argument("memuse_report", nargs='+',
                    help="Memory usage reports (order by step if using more than one)")
parser.add_argument(
    "-b",
    "--blacklist",
    dest="blacklist",
    help="substring of allocations to ignore (maybe be repeated)",
    default=None,
    action='append'
)
args = parser.parse_args()

memuse = OrderedDict()
labels = {}
totalmem = 0
process_use = ""
peak = 0.0

for filename in args.memuse_report:
    sys.stderr.write("## Processing: " + filename + "\n")
    with open(filename) as infile:
        header = np.load(infile)
        process_use = header['metadata'][0][14:].split('\n')[0]
        usedlabels = np.load(infile)
        logdata = np.load(infile)
        print '# {:<18s} {:>30s} {:>9s} {:>9s} {:s}'.format("tic", "label", "allocated", "step", "MB")
        for line in np.nditer(logdata):
            step = int(line['step'])
            allocated = int(line['allocated'])
            size = int(line['size'])
            adr = int(line['ptr'])
            tic = int(line['dtic'])
            label = str(usedlabels[line['index']]['label'])

            #  Skip blacklisted allocations, these can swamp the signal...
            if args.blacklist != None:
                skip = False
                for item in args.blacklist:
                    if item in label:
                        skip = True
                        break
                if skip:
                    continue

            doprint = True
            if allocated == 1:
                #  Allocation.
                totalmem = totalmem + size
                if not adr in memuse:
                    memuse[adr] = [size]
                    labels[adr] = label
                else:
                    memuse[adr].append(size)
            else:
                #  Free, locate allocation.
                if not adr in memuse:
                    #  Unmatched free, complain and skip.
                    #print "### unmatched free: ", label, adr
                    doprint = False
                else:
                    allocs = memuse[adr]
                    totalmem = totalmem - allocs[0]
                    if len(allocs) > 1:
                        memuse[adr] = allocs[1:]
                    else:
                        del memuse[adr]
            if doprint:
                if totalmem > peak:
                    peak = totalmem
                print '{:<20d} {:>30s} {:9d} {:9d} {:.3f}'.format(tic, label, allocated, step, totalmem/(1048576.0))
    sys.stderr.write("## Finished ingestion of: " + filename + "\n")

totals = {}
numactive = {}
for adr in labels:
    #  If any remaining allocations.
    if adr in memuse:
        if labels[adr] in totals:
            totals[labels[adr]] = totals[labels[adr]] + memuse[adr][0]
            numactive[labels[adr]] = numactive[labels[adr]] + 1
        else:
            totals[labels[adr]] = memuse[adr][0]
            numactive[labels[adr]] = 1

print "# Memory use by label:"
print "## ", '{:<30s} {:>16s} {:>16s}'.format("label", "MB", "numactive")
print "## "
total = 0.0
for label in sorted(totals):
    mem = totals[label]/(1048576.0)
    total = total + mem
    print "## ", '{:<30s} {:16.3f} {:16d}'.format(label, mem, numactive[label])
print "## "
print "# Total memory still in use : ", '{:.3f}'.format(total), " (MB)"
print "# Peak memory usage         : ", '{:.3f}'.format(peak/1048576.0), " (MB)"
if process_use != "":
    print "#"
    print "# Memory use by process (all/system):", process_use
sys.exit(0)
