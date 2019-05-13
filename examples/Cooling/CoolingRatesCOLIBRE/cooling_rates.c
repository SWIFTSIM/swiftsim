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

#if defined(COOLING_COLIBRE) && defined(CHEMISTRY_COLIBRE) && \
    defined(GADGET2_SPH)
#include "cooling/COLIBRE/cooling_rates.h"
#include "cooling/COLIBRE/cooling_tables.h"

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
  struct cooling_function_data cooling;
  struct cosmology cosmo;
  const char *parametersFileName = "./cooling_rates.yml";

  /* Initialize CPU frequency, this also starts time. */
  unsigned long long cpufreq = 0;
  clocks_set_cpufreq(cpufreq);

/* Choke on FP-exceptions */
#ifdef HAVE_FE_ENABLE_EXCEPT
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif

  const int npts = 250;  // number of values for the internal energy at which
                         // cooling rate is evaluated

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
        error("invalid option(s) to cooling_rates");
    }

  // Read the parameter file
  if (params == NULL) error("Error allocating memory for the parameter file.");
  message("Reading runtime parameters from file '%s'", parametersFileName);
  parser_read_file(parametersFileName, params);

  // Init units
  units_init_from_params(&us, params, "InternalUnitSystem");
  phys_const_init(&us, params, &internal_const);

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

  // Init cooling
  cooling_init(params, &us, &internal_const, &cooling);
  cooling_print(&cooling);

  // extract mass fractions, calculate table indices and offsets
  float XH = p.chemistry_data.metal_mass_fraction[chemistry_element_H];

  // Calculate contributions from metals to cooling rate
  // open file
  FILE *output_file = fopen("cooling_output.dat", "w");
  if (output_file == NULL) {
    error("Error opening output file!\n");
  }

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
  int red_index, met_index, n_H_index;

  get_index_1d(cooling.Redshifts, colibre_cooling_N_redshifts, cosmo.z,
               &red_index, &d_red);
  get_index_1d(cooling.Metallicity, colibre_cooling_N_metallicity, logZZsol,
               &met_index, &d_met);
  get_index_1d(cooling.nH, colibre_cooling_N_density, log10(inn_h), &n_H_index,
               &d_n_H);

  // Loop over internal energy
  for (int j = 0; j < npts; j++) {

    // Update the particle with the new values
    set_quantities(&p, &xp, &us, &cooling, &cosmo, &internal_const, nh,
                   pow(10.0, 9.0 + j * 9.5 / npts));

    // New internal energy
    u = hydro_get_physical_internal_energy(&p, &xp, &cosmo) *
        cooling.internal_energy_to_cgs;

    // calculate temperature
    const double temperature =
        colibre_convert_u_to_temp(log10(u), cosmo.z, n_H_index, d_n_H,
                                  met_index, d_met, red_index, d_red, &cooling);

    // calculate cooling rates
    const double lambda_net = colibre_cooling_rate(
        log10(u), cosmo.z, nh, pow(10., logZZsol), abundance_ratio, n_H_index,
        d_n_H, met_index, d_met, red_index, d_red, &cooling);

    // Dump...
    fprintf(output_file, "%.5e %.5e\n", exp(M_LN10 * temperature), lambda_net);
  }
  fclose(output_file);
  message("done cooling rates test");

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
