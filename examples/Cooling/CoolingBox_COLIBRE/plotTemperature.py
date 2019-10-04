import h5py 
import numpy as np
import matplotlib
from glob import glob
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import math

# Plot parameters
params = {
    'axes.labelsize': 10,
    'axes.titlesize': 10,
    'font.size': 12,
    'legend.fontsize': 12,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'text.usetex': True,
    'figure.figsize': (5, 5),
    'figure.subplot.left': 0.1,
    'figure.subplot.right': 0.99,
    'figure.subplot.bottom': 0.1,
    'figure.subplot.top': 0.96,
    'figure.subplot.wspace': 0.15,
    'figure.subplot.hspace': 0.12,
    'lines.markersize': 6,
    'lines.linewidth': 3.,
}
plt.rcParams.update(params)

# Physical constants 
m_h_cgs = 1.672621898e-24  # proton mass 

def get_index_1d(table, x): 
    ntable = len(table) 
    denominator = (table[-1] - table[0]) / (ntable - 1) 
    
    if x <= table[0]: 
        i = 0 
        dx = 0.0 
    elif x >= table[ntable - 1]: 
        i = ntable - 2 
        dx = 1.0 
    else: 
        i = int(math.floor((x - table[0]) / denominator)) 
        dx = (x - table[i]) / denominator 

    return i, dx 
    
def interpol_1d(table, i, dx): 
    dx_m = 1.0 - dx 

    output = dx_m * table[i] 
    output += dx * table[i + 1] 

    return output 

def interpol_2d(table, i, j, dx, dy): 
    dx_m = 1.0 - dx 
    dy_m = 1.0 - dy 
    
    output = dx_m * dy_m * table[i, j] 
    output += dx_m * dy * table[i, j + 1] 
    output += dx * dy_m * table[i + 1, j] 
    output += dx * dy * table[i + 1, j + 1]

    return output 

def interpol_3d(table, i, j, k, dx, dy, dz): 
    dx_m = 1.0 - dx 
    dy_m = 1.0 - dy 
    dz_m = 1.0 - dz 

    output = dx_m * dy_m * dz_m * table[i, j, k] 
    output += dx_m * dy_m * dz * table[i, j, k + 1] 
    output += dx_m * dy * dz_m * table[i, j + 1, k] 
    output += dx * dy_m * dz_m * table[i + 1, j, k] 
    output += dx_m * dy * dz * table[i, j + 1, k + 1] 
    output += dx * dy_m * dz * table[i + 1, j, k + 1] 
    output += dx * dy * dz_m * table[i + 1, j + 1, k] 
    output += dx * dy * dz * table[i + 1, j + 1, k + 1] 

    return output 

def temperature_evolution_from_table(cool_table_path, Z_over_Zsol, nH, T_init, XH, time_end): 
    # Read in cooling table. 
    h5file = h5py.File(cool_table_path, "r") 

    U_from_T = np.array(h5file["Tdep/U_from_T"]) 
    cool_rate = np.array(h5file["Tdep/Cooling"]) 
    heat_rate = np.array(h5file["Tdep/Heating"]) 

    # Use the redshift zero rates. 
    redshift_index = 0 

    # Assume solar abundance ratios, so that 
    # we can just use the TotalPrim and 
    # TotalMetals cooling channels. 
    cool_rate_prim = cool_rate[redshift_index, :, :, :, -2] 
    cool_rate_metals = cool_rate[redshift_index, :, :, :, -1] 
    heat_rate_prim = heat_rate[redshift_index, :, :, :, -2] 
    heat_rate_metals = heat_rate[redshift_index, :, :, :, -1] 

    # Get the metallicity and density indices. 
    nH_bins = np.array(h5file["TableBins/DensityBins"]) 
    Z_bins = np.array(h5file["TableBins/MetallicityBins"]) 
    T_bins = np.array(h5file["TableBins/TemperatureBins"]) 

    h5file.close() 

    rho_cgs = nH * m_h_cgs / XH 

    nH_index, dnH = get_index_1d(nH_bins, np.log10(nH)) 
    Z_index, dZ = get_index_1d(Z_bins[1:], np.log10(Z_over_Zsol)) 
    Z_index += 1 

    temperature_array = [T_init] 
    time_array = [0.0] 
    current_time = 0.0 

    # First Step 
    T_index, dT = get_index_1d(T_bins, np.log10(T_init)) 
    u_0 = rho_cgs * (10.0 ** interpol_3d(U_from_T[redshift_index, :, :, :], T_index, Z_index, nH_index, dT, dZ, dnH)) 
    
    # du_dt in erg cm^3 s^-1

    # At start of time-step 
    du_dt_0 = 10.0 ** interpol_3d(heat_rate_prim, T_index, Z_index, nH_index, dT, dZ, dnH) 
    du_dt_0 += 10.0 ** interpol_3d(heat_rate_metals, T_index, Z_index, nH_index, dT, dZ, dnH) 
    du_dt_0 -= 10.0 ** interpol_3d(cool_rate_prim, T_index, Z_index, nH_index, dT, dZ, dnH) 
    du_dt_0 -= 10.0 ** interpol_3d(cool_rate_metals, T_index, Z_index, nH_index, dT, dZ, dnH) 

    # At end of time-step
    if du_dt_0 < 0: 
        # Cooling 
        T_index_1 = T_index 
    else: 
        # Heating 
        T_index_1 = T_index + 1 

    du_dt_1 = 10.0 ** interpol_2d(heat_rate_prim[T_index_1, :, :], Z_index, nH_index, dZ, dnH) 
    du_dt_1 += 10.0 ** interpol_2d(heat_rate_metals[T_index_1, :, :], Z_index, nH_index, dZ, dnH) 
    du_dt_1 -= 10.0 ** interpol_2d(cool_rate_prim[T_index_1, :, :], Z_index, nH_index, dZ, dnH) 
    du_dt_1 -= 10.0 ** interpol_2d(cool_rate_metals[T_index_1, :, :], Z_index, nH_index, dZ, dnH) 

    # Take the average 
    du_dt = (du_dt_0 + du_dt_1) / 2.0 

    # Convert to erg cm^-3 s^-1
    du_dt *= nH ** 2.0 

    if du_dt_0 < 0.0: 
        # Cooling 
        u_1 =  rho_cgs * (10.0 ** interpol_2d(U_from_T[redshift_index, T_index, :, :], Z_index, nH_index, dZ, dnH)) 
        temperature_array.append(10.0 ** T_bins[T_index])
    else: 
        # Heating 
        u_1 = rho_cgs * (10.0 ** interpol_2d(U_from_T[redshift_index, T_index + 1, :, :], Z_index, nH_index, dZ, dnH)) 
        temperature_array.append(10.0 ** T_bins[T_index + 1])

    # Calculate time to cool to 
    # next temperature bin. 
    delta_time = (u_1 - u_0) / du_dt 
    current_time += delta_time 
    time_array.append(current_time)  

    if du_dt < 0.0: 
        # Cool until du_dt becomes positive, 
        # or we reach the minimum temperature 
        # bin, or we reach the end time. 
        while (du_dt < 0.0): 
            if current_time >= time_end: 
                break 

            if (T_index - 1) < 0: 
                break 
              
            u_0 = rho_cgs * (10.0 ** interpol_2d(U_from_T[redshift_index, T_index, :, :], Z_index, nH_index, dZ, dnH)) 
            u_1 = rho_cgs * (10.0 ** interpol_2d(U_from_T[redshift_index, T_index - 1, :, :], Z_index, nH_index, dZ, dnH)) 

            # du_dt at the start of the step 
            du_dt_0 = 10.0 ** interpol_2d(heat_rate_prim[T_index, :, :], Z_index, nH_index, dZ, dnH) 
            du_dt_0 += 10.0 ** interpol_2d(heat_rate_metals[T_index, :, :], Z_index, nH_index, dZ, dnH) 
            du_dt_0 -= 10.0 ** interpol_2d(cool_rate_prim[T_index, :, :], Z_index, nH_index, dZ, dnH) 
            du_dt_0 -= 10.0 ** interpol_2d(cool_rate_metals[T_index, :, :], Z_index, nH_index, dZ, dnH) 

            # du_dt at the end of the step 
            du_dt_1 = 10.0 ** interpol_2d(heat_rate_prim[T_index - 1, :, :], Z_index, nH_index, dZ, dnH) 
            du_dt_1 += 10.0 ** interpol_2d(heat_rate_metals[T_index - 1, :, :], Z_index, nH_index, dZ, dnH) 
            du_dt_1 -= 10.0 ** interpol_2d(cool_rate_prim[T_index - 1, :, :], Z_index, nH_index, dZ, dnH) 
            du_dt_1 -= 10.0 ** interpol_2d(cool_rate_metals[T_index - 1, :, :], Z_index, nH_index, dZ, dnH) 
            
            # Take the average 
            du_dt = (du_dt_0 + du_dt_1) / 2.0 

            # Convert to erg cm^-3 s^-1
            du_dt *= nH ** 2.0 
            
            delta_time = (u_1 - u_0) / du_dt 
            current_time += delta_time 
            time_array.append(current_time) 
            temperature_array.append(10.0 ** T_bins[T_index - 1]) 
            T_index -= 1 
    else: 
        # Heat until du_dt becomes negative, 
        # or we reach the maximum temperature 
        # bin, or we reach the end time. 
        while (du_dt > 0.0): 
            if current_time >= time_end: 
                break 

            if (T_index + 1) >= len(T_bins): 
                break 
                
            u_0 = rho_cgs * (10.0 ** interpol_2d(U_from_T[redshift_index, T_index, :, :], Z_index, nH_index, dZ, dnH)) 
            u_1 = rho_cgs * (10.0 ** interpol_2d(U_from_T[redshift_index, T_index + 1, :, :], Z_index, nH_index, dZ, dnH)) 

            # du_dt at the start of the step 
            du_dt_0 = 10.0 ** interpol_2d(heat_rate_prim[T_index, :, :], Z_index, nH_index, dZ, dnH) 
            du_dt_0 += 10.0 ** interpol_2d(heat_rate_metals[T_index, :, :], Z_index, nH_index, dZ, dnH) 
            du_dt_0 -= 10.0 ** interpol_2d(cool_rate_prim[T_index, :, :], Z_index, nH_index, dZ, dnH) 
            du_dt_0 -= 10.0 ** interpol_2d(cool_rate_metals[T_index, :, :], Z_index, nH_index, dZ, dnH) 

            # du_dt at the end of the step 
            du_dt_1 = 10.0 ** interpol_2d(heat_rate_prim[T_index + 1, :, :], Z_index, nH_index, dZ, dnH) 
            du_dt_1 += 10.0 ** interpol_2d(heat_rate_metals[T_index + 1, :, :], Z_index, nH_index, dZ, dnH) 
            du_dt_1 -= 10.0 ** interpol_2d(cool_rate_prim[T_index + 1, :, :], Z_index, nH_index, dZ, dnH) 
            du_dt_1 -= 10.0 ** interpol_2d(cool_rate_metals[T_index + 1, :, :], Z_index, nH_index, dZ, dnH) 

            # Take the average 
            du_dt = (du_dt_0 + du_dt_1) / 2.0 
            
            # Convert to erg cm^-3 s^-1
            du_dt *= nH ** 2.0 
            
            delta_time = (u_1 - u_0) / du_dt 
            current_time += delta_time 
            time_array.append(current_time) 
            temperature_array.append(10.0 ** T_bins[T_index + 1]) 
            T_index += 1 

    return np.array(time_array), np.array(temperature_array) 

def main(): 
    # Path to COLIBRE cooling table 
    cool_table_path = "/cosma7/data/dp004/dc-ploe1/CoolingTables/2019_most_recent/UV_dust1_CR1_G1_shield1.hdf5" 

    print("Reading snapshot data") 

    # Read the units parameters from the first snapshot
    snap_filename = "output_0000.hdf5"
    f = h5py.File(snap_filename, 'r')
    units = f["InternalCodeUnits"]
    UnitMass_in_cgs = units.attrs["Unit mass in cgs (U_M)"]
    UnitLength_in_cgs = units.attrs["Unit length in cgs (U_L)"]
    UnitTime_in_cgs = units.attrs["Unit time in cgs (U_t)"]

    UnitDensity_in_cgs = UnitMass_in_cgs / (UnitLength_in_cgs ** 3.0) 
    convert_seconds_to_Myr = 1.0 / (60. * 60. * 24. * 365.25 * 1.0e6) 

    # Take simulation properties from the 
    # first snapshot. 
    T = np.array(f["/PartType0/Temperatures"]) 
    rho = np.array(f["/PartType0/Densities"]) 
    metals = np.array(f["/PartType0/ElementMassFractions"]) 
    metallicity = np.array(f["/PartType0/MetalMassFractions"]) 
    XH = metals[0, 0] 
    Z_over_Zsol = metallicity[0] / 0.01337137 
    nH = rho * UnitDensity_in_cgs * XH / m_h_cgs 

    # Take median values 
    ind_sort_T = T.argsort() 
    T_init = T[ind_sort_T][len(T) / 2] 

    ind_sort_nH = nH.argsort() 
    nH_init = nH[ind_sort_nH][len(nH) / 2] 

    f.close() 

    print("T_init = %.4e, nH = %.4e, Z_over_Zsol = %.4e" % (T_init, nH_init, Z_over_Zsol)) 

    # Read snapshots
    files = glob("output_*.hdf5")
    N = len(files)
    T_snap = np.zeros(N)
    time_snap_Myr = np.zeros(N)
    for i in range(N):
        snap = h5py.File(files[i], 'r')
        T = np.array(snap["/PartType0/Temperatures"]) 

        # Take median temperature 
        ind_sort = T.argsort()
        T_sort = T[ind_sort]
        T_median = T_sort[len(T_sort) / 2]
        T_snap[i] = T_median 
        time_snap_Myr[i] = snap["/Header"].attrs["Time"] * UnitTime_in_cgs * convert_seconds_to_Myr
        snap.close()
    
    # Integrate temperature evolution 
    # from cooling tables 
    time_end_Myr = max(time_snap_Myr) 
    time_end_cgs = time_end_Myr / convert_seconds_to_Myr 

    print("Calculating temperature evolution from cooling tables.") 
    time_array, temperature_array = temperature_evolution_from_table(cool_table_path, Z_over_Zsol, nH_init, T_init, XH, time_end_cgs)
    time_array_Myr = time_array * convert_seconds_to_Myr 

    # Scale the y-axis to make it easier 
    # to display the tick labels 
    y_axis_scale = math.floor(np.log10(T_init)) 

    # Plot temperature evolution 
    print("Plotting temperature evolution.") 
    plt.figure()

    # From the simulation 
    plt.plot(time_snap_Myr, T_snap / (10.0 ** y_axis_scale), 'rD', markersize = 3, label = "Simulation")

    # From the cooling table 
    plt.plot(time_array_Myr, temperature_array / (10.0 ** y_axis_scale), 'r-', linewidth = 1.8, label = "Cooling table")

    plt.legend(loc="right", fontsize=8, frameon=False,
               handlelength=3, ncol=1)

    plt.xlim(0.0, time_end_Myr) 
    plt.ylim(0.0, T_init / (10.0 ** y_axis_scale)) 

    y_label = "${\\rm{Temperature~[10^{%d} K]}}$" % (int(y_axis_scale), ) 
    plt.ylabel(y_label, labelpad = 0)
    plt.xlabel("${\\rm{Time~[Myr]}}$", labelpad = 0)

    plt.savefig("temperature_evolution.png", dpi=200)

    plt.close() 

    print("Finished.") 
    
    return 

if __name__ == "__main__": 
    main() 
