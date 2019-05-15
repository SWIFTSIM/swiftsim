import sys
#---------------------------------------------------------------
lognH = float(sys.argv[1])
#---------------------------------------------------------------
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import gridspec
import os

import h5py
import numpy as np
import warnings
warnings.filterwarnings('ignore')

Gyr = 1.e9 * 365.25 * 3600. * 24.


# isochoric
filename = 'isochoric_lognH%.2f.dat'%(lognH)
data_isochoric = np.loadtxt(filename)
data_isochoric[:,0] = data_isochoric[:,0] / Gyr

# isobaric
filename = 'isobaric_lognH%.2f.dat'%(lognH)
data_isobaric = np.loadtxt(filename)
data_isobaric[:,0] = data_isobaric[:,0] / Gyr

# isochoric_extraterm
filename = 'isochoric_extraterm_lognH%.2f.dat'%(lognH)
data_isochoric_extraterm = np.loadtxt(filename)
data_isochoric_extraterm[:,0] = data_isochoric_extraterm[:,0] / Gyr

titlestring = 'Cooling evolution log nH [ccm] = %.2f'%(lognH)

fig = plt.figure()
#fig.set_size_inches(10,6.2,forward=True)
fig.suptitle(titlestring)
#fig.subplots_adjust(left = 0.1, right = 0.95, bottom = 0.15, top = 0.7)
gs = gridspec.GridSpec(1,1,wspace=0, hspace=0)


ax = plt.subplot(gs[0])

ax.set_xlabel('Time [Gyr]')
ax.set_ylabel('log T [K]')
ax.plot(data_isochoric[:,0], data_isochoric[:,3], label = 'Isochoric')
ax.plot(data_isobaric[:,0], data_isobaric[:,3], label = 'Isobaric')
ax.plot(data_isochoric_extraterm[:,0], data_isochoric_extraterm[:,3], label = 'Isochoric extra term')

handles, labels = ax.get_legend_handles_labels()
ax.legend(handles, labels, loc=1)

plt.show()

