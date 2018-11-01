import numpy as np
import matplotlib.pyplot as plt
from astropy import units
import h5py as h5


def read_data():
    R = np.array([])
    z = np.array([])
    x = np.array([])
    y = np.array([])
    for frame in range(0, 2001, 1):
        try:
            sim = h5.File("output_%04d.hdf5" % frame, "r")
        except IOError:
            break

        boxSize = sim["/Header"].attrs["BoxSize"][0]
        pos = sim["/PartType1/Coordinates"][:, :] - boxSize/2.0
        R = np.append(R, np.sqrt(pos[0, 0]**2 + pos[0, 1]**2))
        z = np.append(z, pos[0, 2])
        x = np.append(x, pos[0, 0])
        y = np.append(y, pos[0, 1])
    return (R, z, x, y)





(R, z, x, y) = read_data()

plt.plot(z,y)
plt.savefig('zyplot.png')
plt.close()
plt.plot(x,y)
plt.savefig('xyplot.png')
plt.close()
plt.plot(z,x)
plt.savefig('zxplot.png')
plt.close()
