
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
			struct swift_params *params) {
  dustevo_properties->with_sputtering =
      parser_get_param_int(params, "COLIBREDustEvolution:use_sputtering");
  dustevo_properties->with_SNII_destruction =
      parser_get_param_int(params, "COLIBREDustEvolution:use_SNII_destruction");
  dustevo_properties->with_accretion =
      parser_get_param_int(params, "COLIBREDustEvolution:use_accretion");
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

void dustevo_struct_restore(const struct dustevo_props* dustevo, FILE* stream) {
  restart_read_blocks((void*)dustevo, sizeof(struct dustevo_props), 1, stream,
                      NULL, "dustevo function");
}

void dustevo_struct_dump(const struct dustevo_props* dustevo,
                           FILE* stream) {
  restart_write_blocks((void*)dustevo, sizeof(struct dustevo_props),
                       1, stream, "dustevo", "dustevo function");
}
