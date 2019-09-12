import h5py
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np
from sphviewer.tools import QuickView
from matplotlib.colors import LogNorm
import glob
import re
cmap_gas = mpl.cm.Greys
cmap_HII = mpl.cm.Reds

###################################
filename = 'output_0010.hdf5'
###################################

kpc_in_cm = 3.086e21
Myr_in_s  = 1.e6 * 365.25 * 24. * 3600 

filename_init = 'output_0000.hdf5'
with h5py.File(filename_init, 'r') as f:
        pgas = f['PartType0/Coordinates'][:]
        unit_length_in_cgs = f["/Units"].attrs["Unit length in cgs (U_L)"]
        unit_mass_in_cgs = f["/Units"].attrs["Unit mass in cgs (U_M)"]
        unit_time_in_cgs = f["/Units"].attrs["Unit time in cgs (U_t)"]

xc = np.median(pgas[:,0])
yc = np.median(pgas[:,1])
zc = np.median(pgas[:,2])

r_img = 15.
xmin = -r_img
ymin = -r_img
xmax =  r_img
ymax =  r_img


reg_exp = "output_*.hdf5"
snaplist = sorted(glob.glob(reg_exp))
#######################################################################
for snap in snaplist[:-1]:
    filename = snap
    tmp = re.findall(r'\d+', filename)
    isnap = list(map(int, tmp))[0]

    with h5py.File(filename, 'r') as f:
            time = f['Header'].attrs["Time"]
            pgas = f['PartType0/Coordinates'][:]
            mgas = f['PartType0/Masses'][:]
            HIIr = f['PartType0/HIIregionsEndTime'][:]
            hsml = f['PartType0/SmoothingLengths'][:]
            pstar = f['PartType4/Coordinates'][:]
            birth = f['PartType4/BirthTimes'][:]

    xc = np.median(pgas[:,0])
    yc = np.median(pgas[:,1])
    zc = np.median(pgas[:,2])
    time_in_Myr = time * unit_time_in_cgs / Myr_in_s
    hsml_in_kpc = hsml * unit_length_in_cgs / kpc_in_cm

    indxHII = np.where(HIIr > 0.)[0]
    indxyoung = np.where(birth > 0.)[0]

    pgas_cen = np.zeros_like(pgas)
    pgas_cen[:,0] = (pgas[:,0] - xc) * unit_length_in_cgs / kpc_in_cm
    pgas_cen[:,1] = (pgas[:,1] - yc) * unit_length_in_cgs / kpc_in_cm
    pgas_cen[:,2] = (pgas[:,2] - zc) * unit_length_in_cgs / kpc_in_cm

    pstar_cen = np.zeros_like(pstar)
    pstar_cen[:,0] = (pstar[:,0] - xc) * unit_length_in_cgs / kpc_in_cm
    pstar_cen[:,1] = (pstar[:,1] - yc) * unit_length_in_cgs / kpc_in_cm
    pstar_cen[:,2] = (pstar[:,2] - zc) * unit_length_in_cgs / kpc_in_cm

    ###### Normal gas plot ########################
    qv_gas = QuickView(pgas_cen, hsml = hsml_in_kpc, mass=np.ones(len(pgas_cen)), \
                       plot=False, r='infinity', p=0, t=180, extent =[xmin, xmax, ymin, ymax], 
                       x = 0, y = 0, z = 0)
    img_gas = qv_gas.get_image()
    ext_gas = qv_gas.get_extent()

    ###### HII plot ########################
    if len(indxHII) > 0:
            qv_HII = QuickView(pgas_cen[indxHII], hsml = hsml_in_kpc[indxHII], mass=np.ones(len(pgas_cen[indxHII])), \
                               plot=False, r='infinity', p=0, t=180, extent =[xmin, xmax, ymin, ymax],
                               x = 0, y = 0, z = 0)
            img_HII = qv_HII.get_image()
            ext_HII = qv_HII.get_extent()
    ###############################################

    fig = plt.figure()
    gs = gridspec.GridSpec(1,2)

    ax = plt.subplot(gs[0])
    ax.set_title('Time = %6.2f Myr'%(time_in_Myr)) 
    ax.tick_params(labelleft = False, labelbottom = False, length = 0)
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    im = ax.imshow(img_gas, cmap=cmap_gas, norm = LogNorm(), extent=ext_gas)
    ax.autoscale(False)
    if len(indxHII) > 0:
        ax.scatter(pgas_cen[indxHII,0], pgas_cen[indxHII, 1], color = 'red', s=2)

    ax = plt.subplot(gs[1])
    ax.set_title(r'Imgbox = (%5.1f kpc)$^2$'%(2.*r_img))
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.tick_params(labelleft = False, labelbottom = False, length = 0)
    ax.set(adjustable='box-forced', aspect = 'equal')
    if len(indxHII) > 0:
        im = ax.imshow(img_HII, cmap=cmap_HII, norm = LogNorm(), extent=ext_HII)
        ax.autoscale(False)

    outfile = 'plots/face_on_%4.4i.png'%(isnap)
    fig.savefig(outfile,dpi = 150)
    print('Saved: %s'%(outfile))
    plt.close('all')
