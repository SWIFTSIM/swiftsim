
/* This file's header */
#include "dust.h"

/* Config parameters. */
#include "../config.h"

/* Some standard headers. */
#include <float.h>
#include <hdf5.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>

/* Local includes. */
#include "adiabatic_index.h"
#include "chemistry.h"
#include "cooling.h"
#include "cooling/COLIBRE/cooling_rates.h"
#include "cooling_struct.h"
#include "cooling/COLIBRE/cooling_tables.h"
#include "entropy_floor.h"
#include "error.h"
#include "exp10.h"
#include "hydro.h"
#include "cooling/COLIBRE/interpolate.h"
#include "io_properties.h"
#include "parser.h"
#include "part.h"
#include "physical_constants.h"
#include "space.h"
#include "units.h"

void dustevo_props_init(struct dustevo_props *dustevo_properties,
			struct swift_params *params,
			const struct phys_const *phys_const) {
  message("Initialising dust parameter values");

  dustevo_properties->with_sputtering =
      parser_get_param_int(params, "COLIBREDustEvolution:use_sputtering");

  dustevo_properties->with_SNII_destruction =
      parser_get_param_int(params, "COLIBREDustEvolution:use_SNII_destruction");

  dustevo_properties->with_accretion =
      parser_get_param_int(params, "COLIBREDustEvolution:use_accretion");

  dustevo_properties->with_cooling_on =
      parser_get_param_int(params, "COLIBREDustEvolution:cooling_on");

  dustevo_properties->specific_numSNII =
    7.039463e-3 / phys_const.const_solar_mass;
}

void dustevo_sputter_part(const struct phys_const *phys_const,
                       const struct unit_system *us,
                       const struct cosmology *cosmo,
                       const struct hydro_props *hydro_properties,
                       const struct entropy_floor_properties *floor_props,
                       const struct cooling_function_data *cooling,
		       /* const struct dustevo_props *dustevo_properties, */
                       struct part *restrict p, struct xpart *restrict xp,
			  const float dt, const float dt_therm, const double time){
  float const *metal_fraction =
      chemistry_get_metal_mass_fraction_for_cooling(p);
  const float fGra = metal_fraction[chemistry_element_Gra];
  const float fSil = metal_fraction[chemistry_element_Sil];
  const float fIde = metal_fraction[chemistry_element_Ide];

  if (fGra+fSil+fIde == 0.) return;
  
  /* surface area of a dust grain, assuming 0.1 micron radius and material density of 
     2 g cm^-3 */
  /* const double dmdr_grain_cgs = 4 * M_PI * 1e-10 * 2.; */
  /* const double m_grain_cgs = (4. * M_PI * 1e-15) / 3; */
  /* const double N_grain_cgs = hydro_get_mass(p)*\/m_grain_cgs; */
  /* const double mass_to_cgs = units_cgs_conversion_factor(us, UNIT_CONV_MASS); */
  
  /* grain radius in cgs */ 
  const double grain_radius_cgs = 1e-5;
  /* specific sputtering rate per grain mass fraction  */
  const double dfdr_cgs = 3./grain_radius_cgs;

  /* Convert Hydrogen mass fraction into Hydrogen number density */
  const float XH = metal_fraction[chemistry_element_H];
  const float rho = hydro_get_physical_density(p, cosmo);
  const double n_H = rho * XH / phys_const->const_proton_mass;
  const double n_H_cgs = n_H * cooling->number_density_to_cgs;
  const double dt_cgs = dt * units_cgs_conversion_factor(us, UNIT_CONV_TIME);

  /* surface area of a dust grain, assuming 0.1 micron radius */
  const float temp = cooling_get_temperature(phys_const, hydro_properties, us, cosmo,
  					       cooling, p, xp);
  const double drdt_grain_cgs = -3.2e-18 * n_H_cgs * pow(pow(2e6/temp, 2.5) + 1, -1);
  const double dfdt_cgs = drdt_grain_cgs * dfdr_cgs;
  const double dfrac = dfdt_cgs * dt_cgs;


  /* message("this particle"); */
  /* message("%e", *fGra); */
  /* message("%e", dt_cgs*dfdt_per_grain_dens_cgs*fSil); */
  if (1. + dfrac < 0.) message("%e", 1. + dfrac); 

  p->chemistry_data.metal_mass_fraction[chemistry_element_Gra] = fGra*max(1. + dfrac, 0.);
  p->chemistry_data.metal_mass_fraction[chemistry_element_Sil] = fSil*max(1. + dfrac, 0.);
  p->chemistry_data.metal_mass_fraction[chemistry_element_Ide] = fIde*max(1. + dfrac, 0.);
  
  /* message("%e", temp); */
  /* message("%e", n_H_cgs); */
  /* message("%e", dmdt_grain_cgs); */
  /* float abundance_ratio[colibre_cooling_N_elementtypes]; */
  /* exit(0); */
  /* const float XH = metal_fraction[chemistry_element_H]; */
}

void dustevo_accretion_part(const struct phys_const *phys_const,
                       const struct unit_system *us,
                       const struct cosmology *cosmo,
                       const struct hydro_props *hydro_properties,
                       const struct entropy_floor_properties *floor_props,
                       const struct cooling_function_data *cooling,
		       /* const struct dustevo_props *dustevo_properties, */
                       struct part *restrict p, struct xpart *restrict xp,
			  const float dt, const float dt_therm, const double time){
 
 
  const float X = p->chemistry_data.metal_mass_fraction[chemistry_element_H];
  const float fGra = p->chemistry_data.metal_mass_fraction[chemistry_element_Gra];
  const float fSil = p->chemistry_data.metal_mass_fraction[chemistry_element_Sil];
  const float fIde = p->chemistry_data.metal_mass_fraction[chemistry_element_Ide];

  /* no dust, no accretion */
  if (fGra+fSil+fIde == 0.) return;

  /* grain radius in cgs */ 
  const double grain_radius_cgs = 1e-5;

  /* Convert Hydrogen mass fraction into Hydrogen number density */
  const float rho = hydro_get_physical_density(p, cosmo);
  const float rho_sg = xp->tracers_data.subgrid_temp;

  const double n_H = rho * X / phys_const->const_proton_mass;
  const double n_H_cgs = n_H * cooling->number_density_to_cgs;
  const double n_H_sg_cgs = n_H_cgs * (rho_sg/rho);

  /* get other relevant particle properties */
  const float Z = p->chemistry_data.metal_mass_fraction_total;
  const float temp = cooling_get_temperature(phys_const, hydro_properties, us, cosmo,
					     cooling, p, xp);
  const float temp_sg = xp->tracers_data.subgrid_temp;

  /* accretion timescale multiplicative terms  (Hirashita & Voshchinnikov 2013), 
     implicitly assuming S_acc = 0.3 */
  const float grain_radius_term = grain_radius_cgs / 1e-5;
  const float Z_term  		= 0.0127 / Z; /* Z_sun accessible somewhere? */
  const float n_H_term		= 1.e3 / n_H_cgs;
  const float T_term		= pow(10./temp, 0.5);

  const float mult_terms 	= grain_radius_term * Z_term * n_H_term * T_term;

  /* conversion factor converted ro CGS from yr */
  const float tau_Gra		= 3.132e15 * mult_terms;
  const float tau_Sil		= 1.8126888 * tau_Gra;
  const float tau_Ide		=  tau_Sil; /* this is probably a bad assumption */

  const double dt_cgs = dt * units_cgs_conversion_factor(us, UNIT_CONV_TIME);

  const double growfrac_Gra = dt_cgs / tau_Gra; 
  const double growfrac_Sil = dt_cgs / tau_Sil; 
  const double growfrac_Ide = dt_cgs / tau_Ide; 

  if (growfrac_Gra > 1) message("large growth fraction => %e", growfrac_Gra); 
  if (0 > 1) message("%e %e", temp_sg/temp, n_H_sg_cgs/n_H_cgs); 
  if (0 > 1) message("%e %e", n_H_cgs, temp_sg); 

  /* message("%e", 1./mult_terms);  */

  p->chemistry_data.metal_mass_fraction[chemistry_element_Gra] *= 1. + growfrac_Gra;
  p->chemistry_data.metal_mass_fraction[chemistry_element_Sil] *= 1. + growfrac_Sil;
  p->chemistry_data.metal_mass_fraction[chemistry_element_Ide] *= 1. + growfrac_Ide;
  
}

void dustevo_accretion_part(const struct phys_const *phys_const,
                       const struct unit_system *us,
                       const struct cosmology *cosmo,
                       const struct hydro_props *hydro_properties,
                       const struct entropy_floor_properties *floor_props,
                       const struct cooling_function_data *cooling,
		       /* const struct dustevo_props *dustevo_properties, */
                       struct part *restrict p, struct xpart *restrict xp,
			  const float dt, const float dt_therm, const double time){
 
 
  const float X = p->chemistry_data.metal_mass_fraction[chemistry_element_H];
  const float fGra = p->chemistry_data.metal_mass_fraction[chemistry_element_Gra];
  const float fSil = p->chemistry_data.metal_mass_fraction[chemistry_element_Sil];
  const float fIde = p->chemistry_data.metal_mass_fraction[chemistry_element_Ide];

  /* no dust, no accretion */
  if (fGra+fSil+fIde == 0.) return;

  /* grain radius in cgs */ 
  const double grain_radius_cgs = 1e-5;

  /* Convert Hydrogen mass fraction into Hydrogen number density */
  const float rho = hydro_get_physical_density(p, cosmo);
  const float rho_sg = xp->tracers_data.subgrid_temp;

  const double n_H = rho * X / phys_const->const_proton_mass;
  const double n_H_cgs = n_H * cooling->number_density_to_cgs;
  const double n_H_sg_cgs = n_H_cgs * (rho_sg/rho);

  /* get other relevant particle properties */
  const float Z = p->chemistry_data.metal_mass_fraction_total;
  const float temp = cooling_get_temperature(phys_const, hydro_properties, us, cosmo,
					     cooling, p, xp);
  const float temp_sg = xp->tracers_data.subgrid_temp;

  /* accretion timescale multiplicative terms  (Hirashita & Voshchinnikov 2013), 
     implicitly assuming S_acc = 0.3 */
  const float grain_radius_term = grain_radius_cgs / 1e-5;
  const float Z_term  		= 0.0127 / Z; /* Z_sun accessible somewhere? */
  const float n_H_term		= 1.e3 / n_H_cgs;
  const float T_term		= pow(10./temp, 0.5);

  const float mult_terms 	= grain_radius_term * Z_term * n_H_term * T_term;

  /* conversion factor converted ro CGS from yr */
  const float tau_Gra		= 3.132e15 * mult_terms;
  const float tau_Sil		= 1.8126888 * tau_Gra;
  const float tau_Ide		=  tau_Sil; /* this is probably a bad assumption */

  const double dt_cgs = dt * units_cgs_conversion_factor(us, UNIT_CONV_TIME);

  const double growfrac_Gra = dt_cgs / tau_Gra; 
  const double growfrac_Sil = dt_cgs / tau_Sil; 
  const double growfrac_Ide = dt_cgs / tau_Ide; 

  if (growfrac_Gra > 1) message("large growth fraction => %e", growfrac_Gra); 
  if (0 > 1) message("%e %e", temp_sg/temp, n_H_sg_cgs/n_H_cgs); 
  if (0 > 1) message("%e %e", n_H_cgs, temp_sg); 

  /* message("%e", 1./mult_terms);  */

  p->chemistry_data.metal_mass_fraction[chemistry_element_Gra] *= 1. + growfrac_Gra;
  p->chemistry_data.metal_mass_fraction[chemistry_element_Sil] *= 1. + growfrac_Sil;
  p->chemistry_data.metal_mass_fraction[chemistry_element_Ide] *= 1. + growfrac_Ide;
  
}

void dustevo_struct_restore(const struct dustevo_props* dustevo, FILE* stream) {
  restart_read_blocks((void*)dustevo, sizeof(struct dustevo_props), 1, stream,
                      NULL, "dustevo function");
}

void dustevo_struct_dump(const struct dustevo_props* dustevo,
                           FILE* stream) {
  restart_write_blocks((void*)dustevo, sizeof(struct dustevo_props),
                       1, stream, "dustevo", "dustevo function");
}
