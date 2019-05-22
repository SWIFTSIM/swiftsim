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

Myr = 1.e6 * 365.25 * 3600. * 24.
Gyr = 1.e9 * 365.25 * 3600. * 24.
kB  = 1.38064852e-16
XH = 0.73
mu = 0.6
mu_neutral = 1.22

# outputfilename
outfilename = 'cooling_evolution_lognH%.2f.png'%(lognH)

# explicit isochoric
filename = 'explicit_isochoric_lognH%.2f.dat'%(lognH)
data_exp_isochoric = np.loadtxt(filename)
data_exp_isochoric[:,1] = data_exp_isochoric[:,1] / Myr

# explicit isobaric
filename = 'explicit_isobaric_lognH%.2f.dat'%(lognH)
data_exp_isobaric = np.loadtxt(filename)
data_exp_isobaric[:,1] = data_exp_isobaric[:,1] / Myr

# explicit isobaric density varies
filename = 'explicit_isobaric_lognH%.2f_densvar.dat'%(lognH)
data_exp_isobaric_vary = np.loadtxt(filename)
data_exp_isobaric_vary[:,1] = data_exp_isobaric_vary[:,1] / Myr

#implicit isochoric
filename = 'implicit_isochoric_lognH%.2f.dat'%(lognH)
data_imp_isochoric = np.loadtxt(filename)
data_imp_isochoric[:,1] = data_imp_isochoric[:,1] / Myr

#implicit isobaric
filename = 'implicit_isobaric_lognH%.2f.dat'%(lognH)
data_imp_isobaric = np.loadtxt(filename)
data_imp_isobaric[:,1] = data_imp_isobaric[:,1] / Myr

#implicit isobaric density varies
filename = 'implicit_isobaric_lognH%.2f_densvar.dat'%(lognH)
data_imp_isobaric_vary = np.loadtxt(filename)
data_imp_isobaric_vary[:,1] = data_imp_isobaric_vary[:,1] / Myr

xxmin = np.amin( (data_exp_isochoric[:,1].min(), data_exp_isobaric[:,1].min(), data_exp_isobaric_vary[:,1].min(), \
                  data_imp_isochoric[:,1].min(), data_imp_isobaric[:,1].min(), data_imp_isobaric_vary[:,1].min() ) )

xxmax = np.amax( (data_exp_isochoric[:,1].max(), data_exp_isobaric[:,1].max(), data_exp_isobaric_vary[:,1].max(), \
                  data_imp_isochoric[:,1].max(), data_imp_isobaric[:,1].max(), data_imp_isobaric_vary[:,1].max() ) )

#TconstL_isochor = np.log10(TconstL_isochor)

titlestring = 'Cooling evolution log nH [ccm] = %.2f'%(lognH)

fig = plt.figure()
fig.set_size_inches(10,11,forward=True)
fig.suptitle(titlestring)
#fig.subplots_adjust(left = 0.1, right = 0.95, bottom = 0.15, top = 0.7)
gs = gridspec.GridSpec(3,2,wspace=0.25, hspace=0.35)

#################################################
ax = plt.subplot(gs[0])

ax.set_title('Isochoric: high T')
ax.set_xlabel('Time [Myr]')
ax.set_ylabel('T [K]')
ax.set_xlim(xxmin, xxmax)
ax.set_ylim(bottom = 1.e6, top = 1.1e9)
ax.set_yscale("log", nonposy = 'clip')

ax.plot(data_exp_isochoric[:,1], np.power(10., data_exp_isochoric[:,4]), label = 'Explicit', color = 'black', linestyle = 'solid')
ax.plot(data_imp_isochoric[:,1], np.power(10., data_imp_isochoric[:,4]), label = 'Implicit', color = 'red' , linestyle = 'dotted')

handles, labels = ax.get_legend_handles_labels()
ax.legend(handles, labels, loc=1)

#################################################
ax = plt.subplot(gs[2])

ax.set_title('Isobaric: high T')
ax.set_xlabel('Time [Myr]')
ax.set_ylabel('T [K]')
ax.set_xlim(xxmin, xxmax)
ax.set_ylim(bottom = 1.e6, top = 1.1e9)
ax.set_yscale("log", nonposy = 'clip')

ax.plot(data_exp_isobaric[:,1], np.power(10., data_exp_isobaric[:,4]), label = 'Explicit (n$_H$ const)', color = 'black', linestyle = 'solid')
ax.plot(data_imp_isobaric[:,1], np.power(10., data_imp_isobaric[:,4]), label = 'Implicit (n$_H$ const)', color = 'red' , linestyle = 'dotted')

#################################################
ax = plt.subplot(gs[4])

ax.set_title('Isobaric (dens vary): high T')
ax.set_xlabel('Time [Myr]')
ax.set_ylabel('T [K]')
ax.set_xlim(xxmin, xxmax)
ax.set_ylim(bottom = 1.e6, top = 1.1e9)
ax.set_yscale("log", nonposy = 'clip')

ax.plot(data_exp_isobaric_vary[:,1], np.power(10., data_exp_isobaric_vary[:,4]), label = 'Explicit (n$_H$ vary)', color = 'black', linestyle = 'dashed')
ax.plot(data_imp_isobaric_vary[:,1], np.power(10., data_imp_isobaric_vary[:,4]), label = 'Implicit (n$_H$ vary)', color = 'red' , linestyle = 'dashed')

#################################################
ax = plt.subplot(gs[1])
yymax = 6.

ax.set_title('Isochoric: low T')
ax.set_xlabel('Time [Myr]')
ax.set_ylabel('T [K]')
ax.set_ylim(bottom = 100., top = pow(10., yymax))
ax.set_yscale("log", nonposy = 'clip')

indx1 = np.where(data_exp_isochoric[:,4] <= yymax)
indx2 = np.where(data_imp_isochoric[:,4] <= yymax)

min1 = np.amin( (data_exp_isochoric[indx1[0][0], 1], data_imp_isochoric[indx2[0][0], 1]) )
max1 = np.amax( (data_exp_isochoric[-1, 1]         , data_imp_isochoric[-1, 1]         ) )

ax.set_xlim(0.99*min1, max1)

ax.plot(data_exp_isochoric[:,1], np.power(10., data_exp_isochoric[:,4]), label = 'Explicit', color = 'black', linestyle = 'solid')
ax.plot(data_imp_isochoric[:,1], np.power(10., data_imp_isochoric[:,4]), label = 'Implicit', color = 'red' , linestyle = 'dotted')

#################################################
ax = plt.subplot(gs[3])
yymax = 6.

ax.set_title('Isobaric: low T')
ax.set_xlabel('Time [Myr]')
ax.set_ylabel('T [K]')
ax.set_yscale("log", nonposy = 'clip')
ax.set_ylim(bottom = 100., top = pow(10., yymax))

indx1 = np.where(data_exp_isobaric[:,4] <= yymax)
indx2 = np.where(data_imp_isobaric[:,4] <= yymax)

min1 = np.amin( (data_exp_isobaric[indx1[0][0], 1], data_imp_isobaric[indx2[0][0], 1] ) )
max1 = np.amax( (data_exp_isobaric[-1, 1], data_imp_isobaric[-1, 1] ) ) 

ax.set_xlim(0.99* min1, max1)

ax.plot(data_exp_isobaric[:,1], np.power(10., data_exp_isobaric[:,4]), label = 'Explicit (n$_H$ const)', color = 'black', linestyle = 'solid')
ax.plot(data_imp_isobaric[:,1], np.power(10., data_imp_isobaric[:,4]), label = 'Implicit (n$_H$ const)', color = 'red' , linestyle = 'dotted')

#################################################
ax = plt.subplot(gs[5])
yymax = 6.

ax.set_title('Isobaric (dens vary): low T')
ax.set_xlabel('Time [Myr]')
ax.set_ylabel('T [K]')
ax.set_yscale("log", nonposy = 'clip')
ax.set_ylim(bottom = 100., top = pow(10., yymax))

indx3 = np.where(data_exp_isobaric_vary[:,4] <= yymax)
indx4 = np.where(data_imp_isobaric_vary[:,4] <= yymax)

min1 = np.amin( (data_exp_isobaric_vary[indx3[0][0], 1], data_imp_isobaric_vary[indx4[0][0], 1]) )
max1 = np.amax( (data_exp_isobaric_vary[-1, 1], data_imp_isobaric_vary[-1, 1]) )

ax.set_xlim(0.999* min1, max1)

ax.plot(data_exp_isobaric_vary[:,1], np.power(10., data_exp_isobaric_vary[:,4]), label = 'Explicit (n$_H$ vary)', color = 'black', linestyle = 'dashed')
ax.plot(data_imp_isobaric_vary[:,1], np.power(10., data_imp_isobaric_vary[:,4]), label = 'Implicit (n$_H$ vary)', color = 'red' , linestyle = 'dashed')

print ('Saving image as: %s'%(outfilename))
fig.savefig(outfilename)
plt.show()

