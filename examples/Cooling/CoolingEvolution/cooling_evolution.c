/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (C) 2015 Matthieu Schaller (matthieu.schaller@durham.ac.uk).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
#include "config.h"

/* Some standard headers. */
#include <fenv.h>
#include <unistd.h>

/* Local headers. */
#include "swift.h"

#if defined(COOLING_COLIBRE) && defined(GADGET2_SPH)
#include "cooling/COLIBRE/cooling_rates.h"
#include "cooling/COLIBRE/cooling_tables.h"

//#define ISOCHORIC 0
//#define ISOBARIC 1
//#define ISOCHORIC_EXTRATERM 2

enum {ISOCHORIC, ISOBARIC, ISOCHORIC_EXTRATERM};


/**
 * @brief Assign particle density and entropy corresponding to the
 * hydrogen number density and internal energy specified.
 *
 * @param p Particle data structure
 * @param cooling Cooling function data structure
 * @param cosmo Cosmology data structure
 * @param internal_const Physical constants data structure
 * @param nh Hydrogen number density (cgs units)
 * @param u Internal energy (cgs units)
 */
void set_quantities(struct part *restrict p, struct xpart *restrict xp,
                    const struct unit_system *restrict us,
                    const struct cooling_function_data *restrict cooling,
                    const struct cosmology *restrict cosmo,
                    const struct phys_const *restrict internal_const, float nh,
                    double u) {

  double hydrogen_number_density =
      nh * pow(units_cgs_conversion_factor(us, UNIT_CONV_LENGTH), 3);
  p->rho = hydrogen_number_density * internal_const->const_proton_mass /
           p->chemistry_data.metal_mass_fraction[chemistry_element_H];

  float pressure = (u * cosmo->a * cosmo->a) *
                   cooling->internal_energy_from_cgs * p->rho *
                   (hydro_gamma_minus_one);
  p->entropy = pressure * (pow(p->rho, -hydro_gamma));
  xp->entropy_full = p->entropy;

  for (int j = 0; j < chemistry_element_count; j++) {
    p->chemistry_data.smoothed_metal_mass_fraction[j] =
        p->chemistry_data.metal_mass_fraction[j];
  }
}

/* 
 * icase = 0 ... isochoric 
 * icase = 1 ... isobaric
 * icase = 2 ... isochoric + extra dN/dt term
 */


INLINE double dTdt(int icase, double hstep_old, double temp,  double redshift, double n_H, float ZZsol, 
            const float abundance_ratio[colibre_cooling_N_elementtypes], int n_H_index,
            float d_n_H, int met_index, float d_met, int red_index, float d_red, 
            const struct cooling_function_data *restrict cooling) {

  const double const_boltzmann_k_cgs = 1.38064852e-16;

  float log10_U = colibre_convert_temp_to_u (temp, redshift, n_H_index, d_n_H,
                                  met_index, d_met, red_index, d_red, cooling);

  double dT, ne, ntotal, lambda_net;

  double s_iso;


  if ((icase == ISOCHORIC) || (icase == ISOCHORIC_EXTRATERM)) {
     s_iso = 0.;
  } else if (icase == ISOBARIC) {
     s_iso = 1.; 
  } else { 
     printf("UNKNOWN CASE"); return -1; 
  } 

  ne = colibre_electron_density(log10_U, redshift, n_H, ZZsol, abundance_ratio, n_H_index,
        d_n_H, met_index, d_met, red_index, d_red, cooling);

  // calculate total number of particles per ccm
  // only atoms and electrons at the moment
  ntotal = 0.;
  for (int ii = 0 ; ii < chemistry_element_count; ii++) {
        ntotal += abundance_ratio[ii];
  }
  ntotal = ne + n_H * ntotal;

  // calculate standard net cooling
  lambda_net = colibre_cooling_rate(
      log10_U, redshift, n_H, ZZsol, abundance_ratio, n_H_index,
      d_n_H, met_index, d_met, red_index, d_red, cooling, 0, 0, 0, 0);

  dT = ne * n_H * lambda_net / ( (3./2. + s_iso) * ntotal * const_boltzmann_k_cgs);

 
  if (icase == ISOCHORIC_EXTRATERM) {
    
    double substep = hstep_old;

    double Told = pow(10., temp);

    if (Told + substep * dT < 10.) printf("temp = %.4f, hstep_old = %.4e, Told = %.4e, Tnew = %.4e\n", temp, hstep_old, Told, Told + substep * dT);

    double Tnew = max(Told + substep * dT, 10.);

    double log10_U_new = colibre_convert_temp_to_u ( log10(Tnew), redshift, n_H_index, d_n_H,
                                 met_index, d_met, red_index, d_red, cooling);
    double ne_new = colibre_electron_density(log10_U_new, redshift, n_H, ZZsol, abundance_ratio, n_H_index,
       d_n_H, met_index, d_met, red_index, d_red, cooling);

    double ntotal_new = ntotal - ne + ne_new;

    printf("dT = %.4e, extra term = %.4e\n", dT, Told / ntotal * (ntotal_new - ntotal) / substep);
    dT = dT - Told / ntotal *  (ntotal_new - ntotal) / substep;
    printf("ntotal = %.4e, ntotal_new = %.4e, hstep_old = %.4e\n", ntotal, ntotal_new, substep);
 
   // dT = dT - pow(10., temp) / ntotal * (ntotal_new - ntotal) / hstep_old;
  }

  return dT;
}

INLINE double make_rkstep(int icase, double hstep, double temp,  double redshift, double n_H, float ZZsol, 
            const float abundance_ratio[colibre_cooling_N_elementtypes], int n_H_index,
            float d_n_H, int met_index, float d_met, int red_index, float d_red, 
            const struct cooling_function_data *restrict cooling, double *hstep_out) {

      double a21 = 1./5.;
      double a31 = 3./40., a32 = 9./40.;
      double a41 = 44./45., a42 = -56./15., a43 = 32./9.;
      double a51 = 19372./6561., a52 = -25360./2187., a53 = 64448./6561., a54 = -212./729.;
      double a61 = 9017./3168., a62 = -355./33., a63 = 46732./5247., a64 = 49./176., a65 = -5103./18656.;

      double b1 = 35./384., b2 = 0., b3 = 500./1113., b4 = 125./192., b5 = -2187./6784., b6 = 11./84.;
      double b1s = 5179./57600., b2s = 0., b3s = 7571./16695., b4s = 393./640., b5s = - 92097./339200., b6s = 187./2100.;

      double rtol = 1.e-3, S = 0.95;

      int count = 0;

      double hstep_old, hstep_new;


      hstep_new = hstep;
      do {
          hstep_old = hstep_new;
          printf("hstep_old = %.4e\n", hstep_old);
          // 6-th order Runge-Kutta
          double tempk1 = temp;
	  printf("T(k1) = %.4f\n", tempk1);
          double k1 = hstep_old * dTdt(icase, hstep_old, tempk1,
            redshift, n_H, ZZsol, abundance_ratio, n_H_index,
            d_n_H, met_index, d_met, red_index, d_red, cooling);
    
          double tempk2 = log10(max(pow(10., temp) + a21 * k1, 10.));
	  printf("T(k2) = %.4f\n", tempk2);
          double k2 = hstep_old * dTdt(icase, hstep_old, tempk2,
            redshift, n_H, ZZsol, abundance_ratio, n_H_index,
            d_n_H, met_index, d_met, red_index, d_red, cooling);
    
          double tempk3 = log10(max(pow(10., temp) + a31 * k1 + a32 * k2, 10.));
	  printf("T(k3) = %.4f\n", tempk3);
          double k3 = hstep_old * dTdt(icase, hstep_old, tempk3,
            redshift, n_H, ZZsol, abundance_ratio, n_H_index,
            d_n_H, met_index, d_met, red_index, d_red, cooling);
    
          double tempk4 = log10(max(pow(10., temp) + a41 * k1 + a42 * k2 + a43 * k3, 10.));
	  printf("T(k4) = %.4f\n", tempk4);
          double k4 = hstep_old * dTdt(icase, hstep_old, tempk4,
            redshift, n_H, ZZsol, abundance_ratio, n_H_index,
            d_n_H, met_index, d_met, red_index, d_red, cooling);
    
          double tempk5 = log10(max(pow(10., temp) + a51 * k1 + a52 * k2 + a53 * k3 + a54 * k4, 10.));
	  printf("T(k5) = %.4f\n", tempk5);
          double k5 = hstep_old * dTdt(icase, hstep_old, tempk5,
            redshift, n_H, ZZsol, abundance_ratio, n_H_index,
            d_n_H, met_index, d_met, red_index, d_red, cooling);
    
          double tempk6 = log10(max(pow(10., temp) + a61 * k1 + a62 * k2 + a63 * k3 + a64 * k4 + a65 * k5, 10.));
	  printf("T(k6) = %.4f\n", tempk6);
          double k6 = hstep_old * dTdt(icase, hstep_old, tempk6,
            redshift, n_H, ZZsol, abundance_ratio, n_H_index,
            d_n_H, met_index, d_met, red_index, d_red, cooling);
    
          // lin
          double temperature_new   = pow(10., temp) + b1 *k1 + b2 *k2 + b3 *k3 + b4 *k4 + b5 *k5 + b6 *k6;
          double temperature_news  = pow(10., temp) + b1s*k1 + b2s*k2 + b3s*k3 + b4s*k4 + b5s*k5 + b6s*k6;
      
          printf("temperature_new = %.4e, temperature_news = %.4e\n", temperature_new, temperature_news);

          double error = abs(temperature_new - temperature_news);
          double scale = max(temperature_new, pow(10., temp)) * rtol;
   
          printf("error = %.4e, scale = %.4e\n", error, scale);

          if (error <= scale) { 
            *hstep_out = hstep_old;
            return temperature_new; 
          } else {
            hstep_new = S * hstep_old * pow(1. / error, 1./5.); 
          }
      count += 1;
    } while (count < 1000.);

    printf("RK does not converge: temp = %.4f", temp);
    return -1.;
}


/**
 * @brief Produces contributions to cooling rates for different
 * hydrogen number densities, from different metals,
 * tests 1d and 4d table interpolations produce
 * same results for cooling rate, dlambda/du and temperature.
 */
int main(int argc, char **argv) {
  // Declare relevant structs
  struct swift_params *params = malloc(sizeof(struct swift_params));
  struct unit_system us;
  struct chemistry_global_data chem_data;
  struct part p;
  struct xpart xp;
  struct phys_const internal_const;
  struct hydro_props hydro_properties;
  struct cooling_function_data cooling;
  struct cosmology cosmo;
  const char *parametersFileName = "./cooling_evolution.yml";
  double hstep_out;

  /* Initialize CPU frequency, this also starts time. */
  unsigned long long cpufreq = 0;
  clocks_set_cpufreq(cpufreq);

/* Choke on FP-exceptions */
#ifdef HAVE_FE_ENABLE_EXCEPT
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif

  // Set some default values
  float redshift = 0.0, log_10_nh = -1;

  // Read options
  int param;
  while ((param = getopt(argc, argv, "z:d:")) != -1) switch (param) {
      case 'z':
        // read redshift
        redshift = atof(optarg);
        break;
      case 'd':
        // read log10 of hydrogen number density
        log_10_nh = atof(optarg);
        break;
      case '?':
        if (optopt == 'z')
          printf("Option -%c requires an argument.\n", optopt);
        else
          printf("Unknown option character `\\x%x'.\n", optopt);
        error("invalid option(s) to cooling_evolution");
    }

  // Read the parameter file
  if (params == NULL) error("Error allocating memory for the parameter file.");
  message("Reading runtime parameters from file '%s'", parametersFileName);
  parser_read_file(parametersFileName, params);

  // Init units
  units_init_from_params(&us, params, "InternalUnitSystem");
  phys_const_init(&us, params, &internal_const);

  // Init properties of hydro
  hydro_props_init(&hydro_properties, &internal_const, &us, params);

  // Init chemistry
  chemistry_init(params, &us, &internal_const, &chem_data);
  chemistry_first_init_part(&internal_const, &us, &cosmo, &chem_data, &p, &xp);
  chemistry_print(&chem_data);

  // Init cosmology
  cosmology_init(params, &us, &internal_const, &cosmo);

  // Set redshift and associated quantities
  const float scale_factor = 1.0 / (1.0 + redshift);
  integertime_t ti_current =
      log(scale_factor / cosmo.a_begin) / cosmo.time_base;
  cosmology_update(&cosmo, &internal_const, ti_current);
  message("Redshift is %f", cosmo.z);
  message("Density is [log nH]%f", log_10_nh);

  // Init cooling
  cooling_init(params, &us, &internal_const, &hydro_properties, &cooling);
  cooling_print(&cooling);

  // extract mass fractions, calculate table indices and offsets
  float XH = p.chemistry_data.metal_mass_fraction[chemistry_element_H];

  char isochoric_outputfile[35];
  sprintf(isochoric_outputfile, "isochoric_lognH%.2f.dat", log_10_nh);
  char isobaric_outputfile[35];
  sprintf(isobaric_outputfile, "isobaric_lognH%.2f.dat", log_10_nh);
  char isochoric_extraterm_outputfile[35];
  sprintf(isochoric_extraterm_outputfile, "isochoric_extraterm_lognH%.2f.dat", log_10_nh);

  // Calculate contributions from metals to cooling rate
  // open file
  FILE *ofile_isochoric = fopen(isochoric_outputfile, "w");
  if (ofile_isochoric == NULL) {
    error("Error opening output file!\n");
  }
  FILE *ofile_isobaric = fopen(isobaric_outputfile, "w");
  if (ofile_isobaric == NULL) {
    error("Error opening output file (heat)!\n");
  }
  FILE *ofile_isochoric_extraterm = fopen(isochoric_extraterm_outputfile, "w");
  if (ofile_isochoric_extraterm == NULL) {
    error("Error opening output file (net)!\n");
  }

  FILE *ofile;

  // set hydrogen number density
  const float nh = exp(M_LN10 * log_10_nh);

  /* Initial internal energy */
  double u = 1.0e14;

  // set internal energy to dummy value, will get reset when looping over
  // internal energies
  set_quantities(&p, &xp, &us, &cooling, &cosmo, &internal_const, nh, u);
  float abundance_ratio[chemistry_element_count + 3];
  float logZZsol = abundance_ratio_to_solar(&p, &cooling, abundance_ratio);

  float inn_h = hydro_get_physical_density(&p, &cosmo) * XH /
                internal_const.const_proton_mass *
                cooling.number_density_to_cgs;

  float d_red, d_met, d_n_H;
  int red_index, met_index, n_H_index, icase;

  double temperature;
  double temperature_new;

  double Tmin = 1.e2;

  get_index_1d(cooling.Redshifts, colibre_cooling_N_redshifts, cosmo.z,
               &red_index, &d_red);
  get_index_1d(cooling.Metallicity, colibre_cooling_N_metallicity, logZZsol,
               &met_index, &d_met);
  get_index_1d(cooling.nH, colibre_cooling_N_density, log10(inn_h), &n_H_index,
               &d_n_H);


  for (icase = 0; icase < 3; icase++) {

      if (icase == ISOCHORIC) {
           ofile = fopen(isochoric_outputfile, "w");
      } else if (icase == ISOBARIC) {
           ofile = fopen(isobaric_outputfile, "w");
      } else if (icase == ISOCHORIC_EXTRATERM) {
           ofile = fopen(isochoric_extraterm_outputfile, "w");       
      } 
      if (ofile == NULL) {
           error("Error opening output file, icase = %i\n", icase);
      }
    
      // Update the particle with the new values
      set_quantities(&p, &xp, &us, &cooling, &cosmo, &internal_const, nh,
                       pow(10.0, 18.0));
    
      // New internal energy
      u = hydro_get_physical_internal_energy(&p, &xp, &cosmo) *
            cooling.internal_energy_to_cgs;
    
      // calculate temperature
      // lin
      temperature = pow(10., 
            colibre_convert_u_to_temp(log10(u), cosmo.z, n_H_index, d_n_H,
                                      met_index, d_met, red_index, d_red, &cooling)); 
    
      // lin
      temperature_new = temperature;
      printf("Starting temperature = %.4e\n", temperature);
    
      double time = 0.;
      do {
        u = pow(10., colibre_convert_temp_to_u (log10(temperature), cosmo.z, n_H_index, d_n_H,
                                      met_index, d_met, red_index, d_red, &cooling));
    
        // Update the particle with the new values
        set_quantities(&p, &xp, &us, &cooling, &cosmo, &internal_const, nh, u);
    
        // New internal energy
        u = hydro_get_physical_internal_energy(&p, &xp, &cosmo) *
            cooling.internal_energy_to_cgs;
    
        // calculate temperature
        // log10
        temperature =
            colibre_convert_u_to_temp(log10(u), cosmo.z, n_H_index, d_n_H,
                                      met_index, d_met, red_index, d_red, &cooling);
    
        double hstep = 1.e15;
    
        do {
          temperature_new = make_rkstep(icase, hstep, temperature,
            cosmo.z, nh, pow(10., logZZsol), abundance_ratio, n_H_index,
            d_n_H, met_index, d_met, red_index, d_red, &cooling, &hstep_out);
    
          hstep = hstep / 2.;
    
        } while (temperature_new < 0.);
        printf("Outside of make_rkstep: temperature_new = %.4e\n", temperature_new);
     
        time += hstep_out;
    
    
        fprintf(ofile, "%.4e\t%.4e\t%.4f\t%.4f\n", time, hstep_out, temperature, log10(temperature_new) );
    
        if (temperature_new >= (1. - 1.e-5) * pow(10., temperature)) break;
        // lin
        temperature = temperature_new;
    
      } while( temperature_new > Tmin);
    
    
      fclose(ofile);

  }


  message("done cooling evolution test");

  /* Clean everything */
  cosmology_clean(&cosmo);
  cooling_clean(&cooling);

  free(params);
  return 0;
}

#else

int main(int argc, char **argv) {

  /* Initialize CPU frequency, this also starts time. */
  unsigned long long cpufreq = 0;
  clocks_set_cpufreq(cpufreq);

  message(
      "This test is only defined for the COLIBRE cooling model and Gadget-2 "
      "SPH.");
  return 0;
}
#endif
