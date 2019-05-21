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
#include "cooling/COLIBRE/cooling.c"

enum {ISOCHORIC, ISOBARIC, ISOBARIC_DENSVAR};

INLINE double bisection_iter_isobaric(
    const double u_ini_cgs, const double n_H_cgs, const double redshift,
    int n_H_index, float d_n_H, int met_index, float d_met, int red_index,
    float d_red, double Lambda_He_reion_cgs, double ratefact_cgs, float ZZsol,
    const struct cooling_function_data *restrict cooling,
    const float abundance_ratio[chemistry_element_count + 3], double dt_cgs,
    long long ID) {

  /* Bracketing */
  double u_lower_cgs = u_ini_cgs;
  double u_upper_cgs = u_ini_cgs;

  int icoolcase = 0;

  /*************************************/
  /* Let's get a first guess           */
  /*************************************/

  double LambdaNet_cgs =
      Lambda_He_reion_cgs +
      colibre_cooling_rate(log10(u_ini_cgs), redshift, n_H_cgs, ZZsol,
                           abundance_ratio, n_H_index, d_n_H, met_index, d_met,
                           red_index, d_red, cooling, icoolcase, icoolcase,
                           icoolcase, icoolcase);

  LambdaNet_cgs = 3./5. * LambdaNet_cgs;

  /*************************************/
  /* Let's try to bracket the solution */
  /*************************************/

  if (LambdaNet_cgs < 0) {

    /* we're cooling! */
    u_lower_cgs /= bracket_factor;
    u_upper_cgs *= bracket_factor;

    /* Compute a new rate */
    LambdaNet_cgs =
        Lambda_He_reion_cgs +
        colibre_cooling_rate(log10(u_lower_cgs), redshift, n_H_cgs, ZZsol,
                             abundance_ratio, n_H_index, d_n_H, met_index,
                             d_met, red_index, d_red, cooling, icoolcase,
                             icoolcase, icoolcase, icoolcase);
    LambdaNet_cgs = 3./5. * LambdaNet_cgs;

    int i = 0;
    while (u_lower_cgs - u_ini_cgs - LambdaNet_cgs * ratefact_cgs * dt_cgs >
               0 &&
           i < bisection_max_iterations) {

      u_lower_cgs /= bracket_factor;
      u_upper_cgs /= bracket_factor;

      /* Compute a new rate */
      LambdaNet_cgs =
          Lambda_He_reion_cgs +
          colibre_cooling_rate(log10(u_lower_cgs), redshift, n_H_cgs, ZZsol,
                               abundance_ratio, n_H_index, d_n_H, met_index,
                               d_met, red_index, d_red, cooling, icoolcase,
                               icoolcase, icoolcase, icoolcase);
      LambdaNet_cgs = 3./5. * LambdaNet_cgs;
      i++;
    }
    if (i >= bisection_max_iterations) {
      error(
          "particle %llu exceeded max iterations searching for bounds when "
          "cooling",
          ID);
    }
  } else {

    /* we are heating! */
    u_lower_cgs /= bracket_factor;
    u_upper_cgs *= bracket_factor;

    /* Compute a new rate */
    LambdaNet_cgs =
        Lambda_He_reion_cgs +
        colibre_cooling_rate(log10(u_upper_cgs), redshift, n_H_cgs, ZZsol,
                             abundance_ratio, n_H_index, d_n_H, met_index,
                             d_met, red_index, d_red, cooling, icoolcase,
                             icoolcase, icoolcase, icoolcase);
    LambdaNet_cgs = 3./5. * LambdaNet_cgs;

    int i = 0;
    while (u_upper_cgs - u_ini_cgs - LambdaNet_cgs * ratefact_cgs * dt_cgs <
               0 &&
           i < bisection_max_iterations) {

      u_lower_cgs *= bracket_factor;
      u_upper_cgs *= bracket_factor;

      /* Compute a new rate */
      LambdaNet_cgs =
          Lambda_He_reion_cgs +
          colibre_cooling_rate(log10(u_upper_cgs), redshift, n_H_cgs, ZZsol,
                               abundance_ratio, n_H_index, d_n_H, met_index,
                               d_met, red_index, d_red, cooling, icoolcase,
                               icoolcase, icoolcase, icoolcase);
      LambdaNet_cgs = 3./5. * LambdaNet_cgs;
      i++;
    }

    if (i >= bisection_max_iterations) {
      error(
          "particle %llu exceeded max iterations searching for bounds when "
          "heating",
          ID);
    }
  }

  /********************************************/
  /* We now have an upper and lower bound.    */
  /* Let's iterate by reducing the bracketing */
  /********************************************/
  /* bisection iteration */
  int i = 0;
  double u_next_cgs;

  do {

    /* New guess */
    u_next_cgs = 0.5 * (u_lower_cgs + u_upper_cgs);

    /* New rate */
    LambdaNet_cgs =
        Lambda_He_reion_cgs +
        colibre_cooling_rate(log10(u_next_cgs), redshift, n_H_cgs, ZZsol,
                             abundance_ratio, n_H_index, d_n_H, met_index,
                             d_met, red_index, d_red, cooling, icoolcase,
                             icoolcase, icoolcase, icoolcase);
    LambdaNet_cgs = 3./5. * LambdaNet_cgs;

    /* Where do we go next? */
    if (u_next_cgs - u_ini_cgs - LambdaNet_cgs * ratefact_cgs * dt_cgs > 0.0) {
      u_upper_cgs = u_next_cgs;
    } else {
      u_lower_cgs = u_next_cgs;
    }

    i++;
  } while (fabs(u_upper_cgs - u_lower_cgs) / u_next_cgs > bisection_tolerance &&
           i < bisection_max_iterations);

  if (i >= bisection_max_iterations)
    error("Particle id %llu failed to converge", ID);

  return u_upper_cgs;
}




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



INLINE double dTdt_cooling(float XH, int icase, double log10_T,  double redshift, double n_H, float ZZsol,
            const float abundance_ratio[colibre_cooling_N_elementtypes], int n_H_index,
            float d_n_H, int met_index, float d_met, int red_index, float d_red,
            const struct cooling_function_data *restrict cooling) {
  
  double s_iso, mu, ntotal, lambda_net, dQ, dTdt;
  const double const_boltzmann_k_cgs = 1.38064852e-16;
  const double m_u = 1.66054e-24;


  if (icase == ISOCHORIC) {
     s_iso = 0.;
  } else if ( (icase == ISOBARIC) || (icase == ISOBARIC_DENSVAR) ) {
     s_iso = 1.;
  } else {
     printf("UNKNOWN CASE\n"); return +1.;
  }

  mu = colibre_meanparticlemass_temperature(log10_T, redshift, n_H, ZZsol, n_H_index,
       d_n_H, met_index, d_met, red_index, d_red, cooling);

  // also takes molecules into account
  ntotal = n_H / (XH * mu);

  // calculate standard net cooling
  lambda_net = colibre_cooling_rate_temperature(
      log10_T, redshift, n_H, ZZsol, abundance_ratio, n_H_index,
      d_n_H, met_index, d_met, red_index, d_red, cooling, 0, 0, 0, 0);

  dQ = n_H * n_H * lambda_net / (ntotal * mu * m_u);
  dTdt = m_u * mu / ((3./2. + s_iso) * const_boltzmann_k_cgs) * dQ;

  return dTdt;
}

INLINE void temperature_evolution_implicit(int icase, float inn_h, double redshift, int red_index, float d_red, 
                                  double XH, float logZZsol, int met_index, float d_met, 
				  const float abundance_ratio[colibre_cooling_N_elementtypes],
                                  const struct cooling_function_data *restrict cooling) {

  double temperature_start = 1.e9;
  double u_new_cgs;
  double logustart;
  int T_index, n_H_index, U_index;
  float d_T, d_n_H, d_U;
  const double Lambda_He_reion_cgs = 0.;
  const double dt_cgs = 1.e12;
  const long pid = 0;
  double ratefact_cgs, mu, pres_start;
  const double const_boltzmann_k_cgs = 1.38064852e-16;
  double ZZsol = pow(10., logZZsol);


  FILE *ofile;

  if (icase == ISOCHORIC) {
           char isochoric_outputfile[35];
           sprintf(isochoric_outputfile, "implicit_isochoric_lognH%.2f.dat", log10(inn_h));
           ofile = fopen(isochoric_outputfile, "w");
  } else if (icase == ISOBARIC) {
           char isobaric_outputfile[35];
           sprintf(isobaric_outputfile , "implicit_isobaric_lognH%.2f.dat", log10(inn_h));
           ofile = fopen(isobaric_outputfile, "w");
  } else if (icase == ISOBARIC_DENSVAR) {
           char isobaric_outputfile_densvar[35];
           sprintf(isobaric_outputfile_densvar , "implicit_isobaric_lognH%.2f_densvar.dat", log10(inn_h));
           ofile = fopen(isobaric_outputfile_densvar, "w");
  }
  if (ofile == NULL) {
           error("Error opening output file, icase = %i\n", icase);
  }

  get_index_1d(cooling->Temp, colibre_cooling_N_temperature, log10(temperature_start),
               &T_index, &d_T);

  get_index_1d(cooling->nH, colibre_cooling_N_density, log10(inn_h), &n_H_index,
               &d_n_H);

  ratefact_cgs = inn_h * (XH * cooling->inv_proton_mass_cgs);


  /* Temperature from internal energy */
  logustart = interpolation_4d(
                      cooling->table.U_from_T, red_index, T_index, met_index, n_H_index, d_red,
                      d_T, d_met, d_n_H, colibre_cooling_N_redshifts,
                      colibre_cooling_N_temperature, colibre_cooling_N_metallicity,
                      colibre_cooling_N_density);

  if (icase == ISOBARIC_DENSVAR) {
    mu = colibre_meanparticlemass_temperature(log10(temperature_start), redshift, inn_h, ZZsol,
                        n_H_index, d_n_H, met_index, d_met, red_index, d_red, cooling);

    pres_start = inn_h * const_boltzmann_k_cgs * temperature_start / ( XH * mu );
  }

   
  double n_H = inn_h;
  double u_0_cgs = pow(10.,logustart);
  double log10_T = log10(temperature_start);
  double log10_Tnew;
  double time = 0.;
  double dTdt;

  do {

     get_index_1d(cooling->Temp, colibre_cooling_N_temperature, log10_T,
               &T_index, &d_T);

     // Density varies
     if (icase == ISOBARIC_DENSVAR) {
         get_index_1d(cooling->nH, colibre_cooling_N_density, log10(n_H), &n_H_index,
                &d_n_H);
         mu = colibre_meanparticlemass_temperature(log10_T, redshift, n_H, ZZsol,
                        n_H_index, d_n_H, met_index, d_met, red_index, d_red, cooling);
         n_H = pres_start * XH * mu / ( const_boltzmann_k_cgs * pow(10., log10_T) );
         get_index_1d(cooling->nH, colibre_cooling_N_density, log10(n_H), &n_H_index,
                &d_n_H);
         mu = colibre_meanparticlemass_temperature(log10_T, redshift, n_H, ZZsol,
                        n_H_index, d_n_H, met_index, d_met, red_index, d_red, cooling);
         n_H = pres_start * XH * mu / ( const_boltzmann_k_cgs * pow(10., log10_T) );
         get_index_1d(cooling->nH, colibre_cooling_N_density, log10(n_H), &n_H_index,
                &d_n_H);
         ratefact_cgs = n_H * (XH * cooling->inv_proton_mass_cgs);
     }

     if (icase == ISOCHORIC) {
          u_new_cgs = bisection_iter(u_0_cgs, n_H, redshift, n_H_index, d_n_H,
                               met_index, d_met, red_index, d_red,
                               Lambda_He_reion_cgs, ratefact_cgs, logZZsol,
                               cooling, abundance_ratio, dt_cgs, pid);
     } else { 
          u_new_cgs = bisection_iter_isobaric(u_0_cgs, n_H, redshift, n_H_index, d_n_H,
                               met_index, d_met, red_index, d_red,
                               Lambda_He_reion_cgs, ratefact_cgs, logZZsol,
                               cooling, abundance_ratio, dt_cgs, pid);
     }

     if ( u_new_cgs > u_0_cgs ) break;

     get_index_1d(cooling->Therm, colibre_cooling_N_internalenergy, log10(u_new_cgs),
                 &U_index, &d_U);

     log10_Tnew = interpolation_4d(
                      cooling->table.T_from_U, red_index, U_index, met_index, n_H_index, d_red,
                      d_U, d_met, d_n_H, colibre_cooling_N_redshifts,
                      colibre_cooling_N_internalenergy, colibre_cooling_N_metallicity,
                      colibre_cooling_N_density);

     time = time + dt_cgs;
     dTdt = ( pow(10., log10_Tnew) - pow(10., log10_T) ) / dt_cgs;
     fprintf(ofile, "%i\t%.4e\t%.4e\t%.4e\t%.4f\t%.4f\n",icase, time, dTdt, dt_cgs, log10_Tnew, log10(n_H));
     u_0_cgs = u_new_cgs;
     log10_T = log10_Tnew;

  } while (log10_T > 2.);


  fclose(ofile);

}

INLINE void temperature_evolution_explicit(int icase, float inn_h, double redshift, int red_index, float d_red, 
                                  double XH, float logZZsol, int met_index, float d_met, 
				  const float abundance_ratio[colibre_cooling_N_elementtypes],
                                  const struct cooling_function_data *restrict cooling) {

   const double eps = 0.05;
   const double const_boltzmann_k_cgs = 1.38064852e-16;
   double temperature_start = 1.e9;

   double time = 0.;
   double n_H = inn_h;
   double dTdt;
   double ZZsol = pow(10., logZZsol);
   float log10_T = log10(temperature_start);
   
   int n_H_index;
   float d_n_H;
   double hstep;
   double pres_start, mu;

   FILE *ofile;

   if (icase == ISOCHORIC) {
           char isochoric_outputfile[35];
           sprintf(isochoric_outputfile, "explicit_isochoric_lognH%.2f.dat", log10(inn_h));
           ofile = fopen(isochoric_outputfile, "w");
   } else if (icase == ISOBARIC) {
           char isobaric_outputfile[35];
           sprintf(isobaric_outputfile , "explicit_isobaric_lognH%.2f.dat", log10(inn_h));
           ofile = fopen(isobaric_outputfile, "w");
   } else if (icase == ISOBARIC_DENSVAR) {
           char isobaric_outputfile_densvar[35];
           sprintf(isobaric_outputfile_densvar , "explicit_isobaric_lognH%.2f_densvar.dat", log10(inn_h));
           ofile = fopen(isobaric_outputfile_densvar, "w");
   }
   if (ofile == NULL) {
           error("Error opening output file, icase = %i\n", icase);
   }


   get_index_1d(cooling->nH, colibre_cooling_N_density, log10(n_H), &n_H_index,
                &d_n_H);

   mu = colibre_meanparticlemass_temperature(log10_T, redshift, n_H,ZZsol,
                        n_H_index, d_n_H, met_index, d_met, red_index, d_red, cooling);

   pres_start = n_H * const_boltzmann_k_cgs * temperature_start / ( XH * mu );


   do {
      if (icase == ISOBARIC_DENSVAR) { 
           mu = colibre_meanparticlemass_temperature(log10_T, redshift, n_H,ZZsol,
                        n_H_index, d_n_H, met_index, d_met, red_index, d_red, cooling);
           n_H = pres_start * XH * mu / ( const_boltzmann_k_cgs * pow(10., log10_T) );
           get_index_1d(cooling->nH, colibre_cooling_N_density, log10(n_H), &n_H_index,
                        &d_n_H);
           mu = colibre_meanparticlemass_temperature(log10_T, redshift, n_H,ZZsol,
                        n_H_index, d_n_H, met_index, d_met, red_index, d_red, cooling);
           n_H = pres_start * XH * mu / ( const_boltzmann_k_cgs * pow(10., log10_T) );
           get_index_1d(cooling->nH, colibre_cooling_N_density, log10(n_H), &n_H_index,
                        &d_n_H);
      }

      dTdt = dTdt_cooling(XH, icase, log10_T, redshift, n_H, ZZsol, abundance_ratio, 
                          n_H_index, d_n_H, met_index, d_met, red_index, d_red, cooling);
   
      if (dTdt > 0.) break;

      hstep = fabs(eps * pow(10., log10_T) / dTdt);
      log10_T = log10( pow(10., log10_T) + hstep * dTdt );
      time = time + hstep;
      fprintf(ofile, "%i\t%.4e\t%.4e\t%.4e\t%.4f\t%.4f\n",icase, time, dTdt, hstep, log10_T, log10(n_H));

   } while (log10_T > 2.);

   fclose(ofile);

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
  printf("XH = %.4f\n", XH);

  // set hydrogen number density
  float nh = exp(M_LN10 * log_10_nh);

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

  float d_red, d_met;
  int red_index, met_index;

  get_index_1d(cooling.Redshifts, colibre_cooling_N_redshifts, cosmo.z,
               &red_index, &d_red);
  get_index_1d(cooling.Metallicity, colibre_cooling_N_metallicity, logZZsol,
               &met_index, &d_met);

  printf("nH/nH = %.4f\t, logZZsol = %.4f\t, inn_h = %.4e\n", abundance_ratio[0], logZZsol, inn_h);

  for (int icase = 0; icase <= ISOBARIC_DENSVAR; icase++) { 
      temperature_evolution_explicit(icase, inn_h, cosmo.z, red_index, d_red, XH, logZZsol, met_index, d_met, 
                                     abundance_ratio, &cooling);
  }
  for (int icase = 0; icase <= ISOBARIC_DENSVAR; icase++) { 
      temperature_evolution_implicit(icase, inn_h, cosmo.z, red_index, d_red, XH, logZZsol, met_index, d_met, 
                                     abundance_ratio, &cooling);
  }


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

