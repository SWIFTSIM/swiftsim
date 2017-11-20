#!/usr/bin/env python3
from pyswiftsim import pointer
from pyswiftsim import cooling

import astropy.units

import matplotlib
import matplotlib.pyplot as plt

import numpy as np

# number of particles
N  = 3000
# params file
params = "/home/loikki/swift_test/cooling_box/coolingBox.yml"

mu = 2

# time step in seconds
dt = 1e1 * astropy.units.yr.to('s')

# read params
params = cooling.read_params(params)
#print(params)

# init units
d = cooling.init_units(params)
units = d["units"]
consts = d["constants"]

dt /= units.UnitTime_in_cgs

# init cooling
cooling_data = cooling.pycooling_init(params, units, consts)

# Init variables
# density
rho = np.logspace(-6, 4, N).astype(np.float) # in cm-3
rho *= consts.const_proton_mass * units.UnitLength_in_cgs**3
# temperature
T = np.logspace(1,9, N).astype(np.float) # in K

# create mesh and reshape
rho, T = np.meshgrid(rho, T)
shape = rho.shape
rho = rho.reshape(np.prod(shape))
T = T.reshape(np.prod(shape))

# u = k_b T / (gamma - 1) mu m_p
# internal energy
u = consts.const_boltzmann_k * T / ((cooling._hydro_gamma - 1.) * mu * consts.const_proton_mass)

# compute cooling
# rate = du / dt
rate = cooling.pycooling_rate(units, cooling_data, consts, rho, u, dt)
# cooling time
t_cool = - u / rate

# sound speed
cs = cooling.soundspeed_from_internal_energy(rho, u)

# reshape
# rho, T
rho = rho.reshape(shape) / (units.UnitLength_in_cgs**3 * consts.const_proton_mass)
T = T.reshape(shape)

# cooling
t_cool = t_cool.reshape(shape) * units.UnitTime_in_cgs
cs = cs.reshape(shape) * units.UnitLength_in_cgs / units.UnitTime_in_cgs

L_cool = np.abs(t_cool * cs / astropy.units.kpc.to("cm"))
print(L_cool.min(), L_cool.max())

plt.title("Cooling Length")
plt.contourf(rho, T, L_cool, locator=matplotlib.ticker.LogLocator())
plt.xlabel("Density [$cm^{-3}$]")
plt.ylabel("Temperature [K]")
c = plt.colorbar()
c.set_label("Cooling Length [kpc]")

ax = plt.gca()
ax.set_xscale('log')
ax.set_yscale('log')


plt.show()
