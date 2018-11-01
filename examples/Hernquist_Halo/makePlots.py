from galpy.potential import NFWPotential
from galpy.potential import HernquistPotential
from galpy.orbit import Orbit
import numpy as np
import matplotlib.pyplot as plt
from astropy import units
import h5py as h5

C = 8.0
M_200 = 2.0
M = 2e12*units.solMass
a = 10*units.kpc


def read_data():
    R = np.array([])
    z = np.array([])
    for frame in range(0, 599, 1):
        try:
            sim = h5.File("output_%04d.hdf5" % frame, "r")
        except IOError:
            break

        boxSize = sim["/Header"].attrs["BoxSize"][0]
        pos = sim["/PartType1/Coordinates"][:, :] - boxSize/2.0
        R = np.append(R, np.sqrt(pos[0, 0]**2 + pos[0, 1]**2))
        z = np.append(z, pos[0, 2])
    return (R, z)



def galpy_nfw_orbit():
    # Setting up the potential
    nfw = HernquistPotential(amp=2*M,a=a)   
    nfw.turn_physical_on()

    vxvv = [50.*units.kpc, 0.*units.km/units.s, 50.*units.km/units.s,
            0.*units.kpc, 0.*units.km/units.s]

    # Calculating the orbit
    ts = np.linspace(0., 0.58, 1000)*units.Gyr
    o = Orbit(vxvv=vxvv)
    o.integrate(ts, nfw, method='odeint')

    return o


o = galpy_nfw_orbit()
(R, z) = read_data()

o.plot()
plt.scatter(R, z, s=5, color='black', marker='x')
plt.savefig("comparison.png")
plt.close()
