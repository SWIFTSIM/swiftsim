import numpy as np
import h5py as h5
import matplotlib.pyplot as plt
import sys

#n_snaps = 101

#some constants
OMEGA = 0.3 # Cosmological matter fraction at z = 0
PARSEC_IN_CGS = 3.0856776e18
KM_PER_SEC_IN_CGS = 1.0e5
CONST_G_CGS = 6.672e-8
h = 0.67777 # hubble parameter
gamma = 5./3.
eta = 1.2349
H_0_cgs = 100. * h * KM_PER_SEC_IN_CGS / (1.0e6 * PARSEC_IN_CGS)
seconds_in_year = 3.154e7
#read some header/parameter information from the first snapshot

filename = "Hydrostatic_000.hdf5"
f = h5.File(filename,'r')
params = f["Parameters"]
unit_mass_cgs = float(params.attrs["InternalUnitSystem:UnitMass_in_cgs"])
unit_length_cgs = float(params.attrs["InternalUnitSystem:UnitLength_in_cgs"])
unit_velocity_cgs = float(params.attrs["InternalUnitSystem:UnitVelocity_in_cgs"])
unit_time_cgs = unit_length_cgs / unit_velocity_cgs
v_c = float(params.attrs["SoftenedIsothermalPotential:vrot"])
v_c_cgs = v_c * unit_velocity_cgs
header = f["Header"]
N = header.attrs["NumPart_Total"][0]
box_centre = np.array(header.attrs["BoxSize"])

#calculate r_vir and M_vir from v_c
r_vir_cgs = v_c_cgs / (10. * H_0_cgs * np.sqrt(OMEGA))
M_vir_cgs = r_vir_cgs * v_c_cgs**2 / CONST_G_CGS

#read the stats output

stats_array = np.genfromtxt("energy.txt",skip_header = 1)
mass_array = stats_array[:,1]
potential_energy_array = stats_array[:,7]/mass_array
internal_energy_array = stats_array[:,4]/mass_array
kinetic_energy_array = stats_array[:,3]/mass_array
radiated_energy_array = stats_array[:,8]/mass_array
time_array = stats_array[:,0]

conversion_factor = unit_velocity_cgs**2 /  v_c_cgs**2
print conversion_factor


# #put into cgs units
# potential_energy_array_cgs = potential_energy_array * unit_velocity_cgs**2 * unit_mass_cgs**2
# kinetic_energy_array_cgs = kinetic_energy_array * unit_velocity_cgs**2 * unit_mass_cgs**2
# internal_energy_array_cgs = internal_energy_array * unit_velocity_cgs**2 * unit_mass_cgs**2
# radiated_energy_array_cgs = radiated_energy_array * unit_velocity_cgs**2 * unit_mass_cgs**2

# put into dimensionless units
pe = potential_energy_array * conversion_factor
ke = kinetic_energy_array * conversion_factor
ie = internal_energy_array * conversion_factor
re = radiated_energy_array * conversion_factor
te = pe + ke + ie + re

#scale relative to initial total energy
# pe /= te[1]
# ke /= te[1]
# ie /= te[1]
# re /= te[1]
# te_change = (te - te[1])/te[1]
# te /= te[1]

#now convert time
time_array_cgs = time_array * unit_time_cgs
dyn_time_cgs = r_vir_cgs / v_c_cgs
t = time_array_cgs / dyn_time_cgs
dyn_time_years = dyn_time_cgs / seconds_in_year
# for i in range(n_snaps):

#     filename = "CoolingHalo_%03d.hdf5" %i
#     f = h5.File(filename,'r')
#     coords_dset = f["PartType0/Coordinates"]
#     coords = np.array(coords_dset)
# #translate coords by centre of box
#     header = f["Header"]
#     snap_time = header.attrs["Time"]
#     snap_time_cgs = snap_time * unit_time_cgs
#     time_array_cgs = np.append(time_array_cgs,snap_time_cgs)
#     coords[:,0] -= box_centre[0]/2.
#     coords[:,1] -= box_centre[1]/2.
#     coords[:,2] -= box_centre[2]/2.
#     radius = np.sqrt(coords[:,0]**2 + coords[:,1]**2 + coords[:,2]**2)
#     radius_cgs = radius*unit_length_cgs
#     radius_over_virial_radius = radius_cgs / r_vir_cgs

#     r = radius_over_virial_radius
#     total_potential_energy = np.sum(v_c**2*np.log(r))
#     potential_energy_array = np.append(potential_energy_array,total_potential_energy)

#     vels_dset = f["PartType0/Velocities"]
#     vels = np.array(vels_dset)
#     speed_squared = vels[:,0]**2 + vels[:,1]**2 + vels[:,2]**2
#     total_kinetic_energy = 0.5 * np.sum(speed_squared)
#     kinetic_energy_array = np.append(kinetic_energy_array,total_kinetic_energy)

#     u_dset = f["PartType0/InternalEnergy"]
#     u = np.array(u_dset)
#     total_internal_energy = np.sum(u)
#     internal_energy_array = np.append(internal_energy_array,total_internal_energy)

# #get the radiated energy

# energy_array = np.genfromtxt("energy.txt",skip_header = 1)
# #rad_energy_time = energy_array[:,0]
# #rad_energy_time_cgs = rad_energy_time * unit_time_cgs
# rad_energy_array = energy_array[:,6]

# #only use every 10th term in the rad_energy_array
# rad_energy_array = rad_energy_array[0::10]

# #put energies in units of v_c^2 and rescale by number of particles

# pe = potential_energy_array / (N*v_c**2)
# ke = kinetic_energy_array / (N*v_c**2)
# ie = internal_energy_array / (N*v_c**2)
# re = rad_energy_array / (N*v_c**2)
# te = pe + ke + ie #+ re

print te


plt.plot(t,ke,label = "Kinetic Energy")
plt.plot(t,pe,label = "Potential Energy")
plt.plot(t,ie,label = "Internal Energy")
plt.plot(t,re,label = "Radiated Energy")
plt.plot(t,te,label = "Total Energy")
plt.legend(loc = "center left")
plt.xlabel(r"$t / t_{dyn}$")
#plt.ylabel(r"$E / v_c^2$")
plt.title(r"$%d \, \, \mathrm{particles} \,,\, v_c = %.1f \, \mathrm{km / s} , \, \, \mathrm{Dynamical\,time} = %e \, \mathrm{years}$" %(N,v_c,dyn_time_years))
#plt.ylabel("Energy per unit mass, relaltive to initial energy")
#plt.ylabel("Fractional change in total energy")
plt.ylabel(r"$\mathrm{Energy\,per\,unit\,mass} / v_c^2$") 
#plt.ylim((-4,2))
#plot_filename = "density_profile_%03d.png" %i
plt.show()

