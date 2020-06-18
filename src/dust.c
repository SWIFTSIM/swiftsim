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
//#include "cooling/COLIBRE/cooling_rates.h"
#include "cooling_struct.h"
//#include "cooling/COLIBRE/cooling_tables.h"
#include "entropy_floor.h"
#include "error.h"
#include "exp10.h"
#include "hydro.h"
//#include "cooling/COLIBRE/interpolate.h"
#include "io_properties.h"
#include "parser.h"
#include "part.h"
#include "physical_constants.h"
#include "space.h"
#include "units.h"

void dustevo_props_init(struct dustevo_props *dustevo_properties,
			struct swift_params *params,
			const struct phys_const *phys_const,
			const struct unit_system *us) {
  message("Initialising dust parameter values");

  dustevo_properties->with_dust_evolution =
      parser_get_param_int(params, "COLIBREDustEvolution:evolve_dust");

  dustevo_properties->with_sputtering =
      parser_get_param_int(params, "COLIBREDustEvolution:use_sputtering");

  dustevo_properties->with_SNII_destruction =
      parser_get_param_int(params, "COLIBREDustEvolution:use_SNII_destruction");

  dustevo_properties->with_accretion =
      parser_get_param_int(params, "COLIBREDustEvolution:use_accretion");

  dustevo_properties->with_subgrid_props =
      parser_get_param_int(params, "COLIBREDustEvolution:use_subgrid_properties");

  dustevo_properties->with_cooling_on =
      parser_get_param_int(params, "COLIBREDustEvolution:cooling_on");

  dustevo_properties->model_type =
      parser_get_param_int(params, "COLIBREDustEvolution:model_type");

  dustevo_properties->pair_to_cooling =
      parser_get_param_int(params, "COLIBREDustEvolution:pair_to_cooling");

  dustevo_properties->specific_numSNII =
    7.039463e-3 / phys_const->const_solar_mass;

  /* dustevo_properties->number_density_to_cgs =  */
  /*   units_cgs_conversion_factor(us, UNIT_CONV_NUMBER_DENSITY); */
  // message("%e", dustevo_properties->number_density_to_cgs);
  // ASK MATTHIEU: why does this return 0 when initialised here??



  dustevo_properties->refractory_idx[0] = chemistry_element_C;
  dustevo_properties->refractory_idx[1] = chemistry_element_O;
  dustevo_properties->refractory_idx[2] = chemistry_element_Mg;
  dustevo_properties->refractory_idx[3] = chemistry_element_Si;
  dustevo_properties->refractory_idx[4] = chemistry_element_Fe;

  dustevo_properties->dust_idx[0] = chemistry_element_Gra;
  dustevo_properties->dust_idx[1] = chemistry_element_Ide;
  dustevo_properties->dust_idx[2] = chemistry_element_Mgd;
  dustevo_properties->dust_idx[3] = chemistry_element_Sil;
  dustevo_properties->dust_idx[4] = chemistry_element_Fed;

  int j;
  /* initialize elements of array n to 0 */         
  for ( j = 0; j < chemistry_element_count; j++ ) {
    dustevo_properties->comp_Gra[j] = 0.;
    dustevo_properties->comp_Sil[j] = 0.;
    dustevo_properties->comp_Ide[j] = 0.;
    dustevo_properties->comp_Mgd[j] = 0.;
    dustevo_properties->comp_Fed[j] = 0.;
  }

  if (dustevo_properties->model_type == 0){
    /* Assume pure graphite grains */
    dustevo_properties->comp_Gra[chemistry_element_C] = -1.;
    dustevo_properties->comp_Gra[chemistry_element_Gra] = 1.;

    /* Assume olivine grains with 1:1 ratio of Mg to Fe grains */
    dustevo_properties->comp_Sil[chemistry_element_O] = -0.3715755;
    dustevo_properties->comp_Sil[chemistry_element_Mg] = -0.1411169;
    dustevo_properties->comp_Sil[chemistry_element_Si] = -0.1630668;
    dustevo_properties->comp_Sil[chemistry_element_Fe] = -0.3242408;
    dustevo_properties->comp_Sil[chemistry_element_Sil] = 1.;

  }

  if (dustevo_properties->model_type == 1){
    /* Carbon content */
    dustevo_properties->comp_Gra[chemistry_element_C] = -1.;
    dustevo_properties->comp_Gra[chemistry_element_Gra] = 1.;
    /* Silicon content */
    dustevo_properties->comp_Sil[chemistry_element_Si] = -1.;
    dustevo_properties->comp_Sil[chemistry_element_Sil] = 1.;
    /* Oxygen content */
    dustevo_properties->comp_Ide[chemistry_element_O] = -1.;
    dustevo_properties->comp_Ide[chemistry_element_Ide] = 1.;
    /* Magnesium content */
    dustevo_properties->comp_Mgd[chemistry_element_Mg] = -1.;
    dustevo_properties->comp_Mgd[chemistry_element_Mgd] = 1.;
    /* Iron content */
    dustevo_properties->comp_Fed[chemistry_element_Fe] = -1.;
    dustevo_properties->comp_Fed[chemistry_element_Fed] = 1.;

  }
  
  /* if (dustevo_properties->model_type == 2){ */
  /*   for () */
  /* } */

}

void evolve_dust_part(const struct phys_const *phys_const,
		      const struct unit_system *us,
		      const struct cosmology *cosmo,
		      const struct hydro_props *hydro_properties,
		      const struct entropy_floor_properties *floor_props,
		      const struct cooling_function_data *cooling,
		      const struct dustevo_props *dp,
		      struct part *restrict p, struct xpart *restrict xp,
		      const float dt, const float dt_therm, const double time){

  const float X = p->chemistry_data.metal_mass_fraction[chemistry_element_H];
  const float Z = p->chemistry_data.metal_mass_fraction_total;
  const float fGra = p->chemistry_data.metal_mass_fraction[chemistry_element_Gra];
  const float fSil = p->chemistry_data.metal_mass_fraction[chemistry_element_Sil];
  const float fIde = p->chemistry_data.metal_mass_fraction[chemistry_element_Ide];

  /* no dust, no accretion */
  if (fGra+fSil+fIde <= 0.) return;

  /* set temperature and density with subgrid or hydro */
  float rho, temp; 

  if (dp->with_subgrid_props) {
    rho = p->cooling_data.subgrid_dens; //xp->tracers_data.subgrid_dens;
    temp = p->cooling_data.subgrid_temp; //xp->tracers_data.subgrid_temp;
  }
  else {
    rho = hydro_get_physical_density(p, cosmo);
    temp = cooling_get_temperature(phys_const, hydro_properties, us, cosmo,
  					       cooling, p, xp);
  }

  /* --------------- Common Parameters --------------- */
  /* grain radius in cgs */ 
  const double grain_radius_cgs = 1e-5;

  /* convert Hydrogen mass fraction into Hydrogen number density */
  const double n_H = rho * X / phys_const->const_proton_mass;


  /* convert some properties to cgs  */
  const double n_H_cgs = n_H * units_cgs_conversion_factor(us, UNIT_CONV_NUMBER_DENSITY);
  const double dt_cgs = dt * units_cgs_conversion_factor(us, UNIT_CONV_TIME);

  /*------------------------- Sputtering -------------------------*/

  double deltf_sput;
  if (dp->with_sputtering){
    /* compute fractional change in grain species mass this timestep */
    const double dfdr_cgs = 3./grain_radius_cgs;
    const double drdt_grain_cgs = -3.2e-18 * n_H_cgs * pow(pow(2e6/temp, 2.5) + 1, -1);
    const double dfdt_cgs = drdt_grain_cgs * dfdr_cgs;
    deltf_sput = dfdt_cgs * dt_cgs;
  }
  else deltf_sput = 0.;

  /*------------------------- SNII shocks -------------------------*/

  double deltf_SNII;

  if (dp->with_SNII_destruction && (xp->sf_data.SFR > 0.)){
    const float eps_eff = 0.1;
    const double unit_mass_cgs = units_cgs_conversion_factor(us, UNIT_CONV_MASS);
    const double unit_time_cgs = units_cgs_conversion_factor(us, UNIT_CONV_TIME);
    /* SNII rate assuming constant SFR, in cgs units */
    const float rate_SNII = xp->sf_data.SFR * dp->specific_numSNII / unit_time_cgs;
  
    /* Yamasawa et. al (2011): the cgs mass of material swept up in SNII shocks */
    const float m_swept = 3.052e+36 * pow(n_H_cgs, -0.202)* pow((Z/0.0127) + 0.039, -0.298);

    /* destruction timescale of dust in cgs */  
    const float tau_dest = p->mass * unit_mass_cgs / (eps_eff*m_swept * rate_SNII);

    /* fractional change in dust due to dust destruction */
    deltf_SNII = -dt_cgs/tau_dest; 
  }
  else {
    deltf_SNII = 0.;
  }

  /* ------------------------- Accretion ------------------------- */

  float deltf_acc_Gra, deltf_acc_Sil, deltf_acc_Ide;

  if (dp->with_accretion){
    /* accretion timescale multiplicative terms  (Hirashita & Voshchinnikov 2013), 
       implicitly assuming S_acc = 0.3 */
    const float grain_radius_term = grain_radius_cgs / 1e-5;
    /* const float Z_term  		= 0.0127 / Z; /\* Z_sun accessible somewhere? *\/ */
    const float n_H_term		= 1.e3 / n_H_cgs;
    const float T_term		= pow(10./temp, 0.5);

    const float mult_terms 	= grain_radius_term * /*Z_term */ n_H_term * T_term;

    /* abundance terms, depending on the key element for each grain species (normalised to solar ratio from Wiersma+09)*/
    /* !!! check consistency with yield factors, better to use an array of sola abundances defined in header? */
    float abundance_ratio_Gra =  2.07e-3/p->chemistry_data.metal_mass_fraction[chemistry_element_C];
    float abundance_ratio_Ide =  1.1e-3/p->chemistry_data.metal_mass_fraction[chemistry_element_Fe];
    float abundance_ratio_Sil =  max(5.91e-4/p->chemistry_data.metal_mass_fraction[chemistry_element_Mg],
					   6.83e-4/p->chemistry_data.metal_mass_fraction[chemistry_element_Si]);
    abundance_ratio_Sil       =  max(1.1e-3/p->chemistry_data.metal_mass_fraction[chemistry_element_Fe],
				     abundance_ratio_Sil);

    /* conversion factor converted to CGS from yr */
    const float tau_Gra		= 3.132e15 * mult_terms * abundance_ratio_Gra;
    const float tau_Sil		= 5.677e15 * mult_terms * abundance_ratio_Sil;
    const float tau_Ide		= 3.132e15 * mult_terms * abundance_ratio_Ide; /* this is probably a bad assumption */

    deltf_acc_Gra = dt_cgs / tau_Gra; 
    deltf_acc_Sil = dt_cgs / tau_Sil; 
    deltf_acc_Ide = dt_cgs / tau_Ide; 
    /* message("deltf_acc_Gra: %e", deltf_acc_Gra); */
  }
  else {
    deltf_acc_Gra = 0.; 
    deltf_acc_Sil = 0.; 
    deltf_acc_Ide = 0.; 
  }


  /* message("%e", max(1. + deltf_sput + deltf_SNII, 0.)); */
  /* compute total fractional change in mass for each grain species */
  double delta_mfrac_Gra = p->chemistry_data.metal_mass_fraction[chemistry_element_Gra] *
    (exp(deltf_sput + deltf_SNII + deltf_acc_Gra) - 1);
  double delta_mfrac_Sil = p->chemistry_data.metal_mass_fraction[chemistry_element_Sil] *
    (exp(deltf_sput + deltf_SNII + deltf_acc_Sil) - 1);
  double delta_mfrac_Ide = p->chemistry_data.metal_mass_fraction[chemistry_element_Ide] *
    (exp(deltf_sput + deltf_SNII + deltf_acc_Ide) - 1);
  
  /* message("%e %e %e", deltf_sput, deltf_SNII, deltf_acc_Gra); */
  /* message("%e %e", exp(deltf_sput+deltf_SNII+deltf_acc_Gra)-1, deltf_sput+deltf_SNII+deltf_acc_Gra); */
  /* if ((1-dp->with_accretion) && (delta_mfrac_Gra )){ */
  /*   message("%e %e", deltf_sput, deltf_SNII); */
  /* } */
  /* message("%e", (1-exp(deltf_sput + deltf_SNII + deltf_acc_Ide))); */
  /* ensure mass growth in dust does not exceed that of current mass in compositional element */
  delta_mfrac_Gra = min(delta_mfrac_Gra,
  			p->chemistry_data.metal_mass_fraction[chemistry_element_C]);
  delta_mfrac_Ide = min(delta_mfrac_Ide,
  			p->chemistry_data.metal_mass_fraction[chemistry_element_Fe]);
  delta_mfrac_Sil = min(delta_mfrac_Sil,
  			p->chemistry_data.metal_mass_fraction[chemistry_element_Mg]);
  delta_mfrac_Sil = min(delta_mfrac_Sil,
  			p->chemistry_data.metal_mass_fraction[chemistry_element_Si]);
  delta_mfrac_Sil = min(delta_mfrac_Sil,
  			p->chemistry_data.metal_mass_fraction[chemistry_element_Fe]);

  /* apply mass fraction changes to elements and grains */
  int i;
  for ( i = 0; i < chemistry_element_count; i++ ) {
    p->chemistry_data.metal_mass_fraction[i] += dp->comp_Gra[i]*delta_mfrac_Gra;   
    p->chemistry_data.metal_mass_fraction[i] += dp->comp_Sil[i]*delta_mfrac_Sil;   
    p->chemistry_data.metal_mass_fraction[i] += dp->comp_Ide[i]*delta_mfrac_Ide;    

    /* if (p->chemistry_data.metal_mass_fraction[i] < 0) { */
    if (0) {
      message("%d: %e ", i, 
	      p->chemistry_data.metal_mass_fraction[i]);
      message("Gra: %e, Sil %e, Ide %e", delta_mfrac_Gra,
	      delta_mfrac_Sil, delta_mfrac_Ide);
    /* message("factors: %e %e %e ", deltf_sput, deltf_SNII, */
    /* 	    deltf_acc_Gra+deltf_acc_Sil+deltf_acc_Ide); */   
      /* exit(0); */
    }

    
  }  
  /* apply mass fraction changes to total metallicity */
  p->chemistry_data.metal_mass_fraction_total -= (delta_mfrac_Gra + delta_mfrac_Sil + delta_mfrac_Ide);
  
}

void evolve_dust_part_m16(const struct phys_const *phys_const,
		      const struct unit_system *us,
		      const struct cosmology *cosmo,
		      const struct hydro_props *hydro_properties,
		      const struct entropy_floor_properties *floor_props,
		      const struct cooling_function_data *cooling,
		      const struct dustevo_props *dp,
		      struct part *restrict p, struct xpart *restrict xp,
		      const float dt, const float dt_therm, const double time){

  const float X = p->chemistry_data.metal_mass_fraction[chemistry_element_H];
  const float Z = p->chemistry_data.metal_mass_fraction_total;
  const float fGra = p->chemistry_data.metal_mass_fraction[chemistry_element_Gra];
  const float fSil = p->chemistry_data.metal_mass_fraction[chemistry_element_Sil];
  const float fIde = p->chemistry_data.metal_mass_fraction[chemistry_element_Ide];
  const float fMgd = p->chemistry_data.metal_mass_fraction[chemistry_element_Mgd];
  const float fFed = p->chemistry_data.metal_mass_fraction[chemistry_element_Fed];


  /* no dust, no accretion */
  if (fGra+fSil+fIde+fMgd+fFed <= 0.) return;

  /* set temperature and density with subgrid or hydro */
  float rho, temp; 

  if (dp->with_subgrid_props) {
    rho = p->cooling_data.subgrid_dens; // xp->tracers_data.subgrid_dens;
    temp = p->cooling_data.subgrid_temp; //xp->tracers_data.subgrid_temp;
  }
  else {
    rho = hydro_get_physical_density(p, cosmo);
    temp = cooling_get_temperature(phys_const, hydro_properties, us, cosmo,
  					       cooling, p, xp);
  }

  /* --------------- Common Parameters --------------- */
  /* grain radius in cgs */ 
  const double grain_radius_cgs = 1e-5;

  /* convert Hydrogen mass fraction into Hydrogen number density */
  const double n_H = rho * X / phys_const->const_proton_mass;

  /* convert some properties to cgs  */
  const double n_H_cgs = n_H * units_cgs_conversion_factor(us, UNIT_CONV_NUMBER_DENSITY);
  const double dt_cgs = dt * units_cgs_conversion_factor(us, UNIT_CONV_TIME);

  /*------------------------- Sputtering -------------------------*/

  double deltf_sput;
  if (dp->with_sputtering){
    /* compute fractional change in grain species mass this timestep */
    const double dfdr_cgs = 3./grain_radius_cgs;
    const double drdt_grain_cgs = -3.2e-18 * n_H_cgs * pow(pow(2e6/temp, 2.5) + 1, -1);
    const double dfdt_cgs = drdt_grain_cgs * dfdr_cgs;
    deltf_sput = dfdt_cgs * dt_cgs;
  }
  else deltf_sput = 0.;

  /*------------------------- SNII shocks -------------------------*/

  double deltf_SNII;

  if (dp->with_SNII_destruction && (xp->sf_data.SFR > 0.)){
    const float eps_eff = 0.1;
    const double unit_mass_cgs = units_cgs_conversion_factor(us, UNIT_CONV_MASS);
    const double unit_time_cgs = units_cgs_conversion_factor(us, UNIT_CONV_TIME);
    /* SNII rate assuming constant SFR, in cgs units */
    const float rate_SNII = xp->sf_data.SFR * dp->specific_numSNII / unit_time_cgs;
  
    /* Yamasawa et. al (2011): the cgs mass of material swept up in SNII shocks */
    const float m_swept = 3.052e+36 * pow(n_H_cgs, -0.202)* pow((Z/0.0127) + 0.039, -0.289);

    /* destruction timescale of dust in cgs */  
    const float tau_dest = p->mass * unit_mass_cgs / (eps_eff*m_swept * rate_SNII);

    /* fractional change in dust due to dust destruction */
    deltf_SNII = -dt_cgs/tau_dest; 
  }
  else {
    deltf_SNII = 0.;
  }

  /* ------------------------- Accretion ------------------------- */

  float deltf_acc_Gra, deltf_acc_Sil, deltf_acc_Ide, deltf_acc_Mgd, deltf_acc_Fed;

  if (dp->with_accretion){
    /* accretion timescale multiplicative terms  (Hirashita & Voshchinnikov 2013), 
       implicitly assuming S_acc = 0.3 */
    const float grain_radius_term = grain_radius_cgs / 1e-5;
    /* const float Z_term  		= 0.0127 / Z; /\* Z_sun accessible somewhere? *\/ */
    const float n_H_term		= 1.e3 / n_H_cgs;
    const float T_term		= pow(10./temp, 0.5);

    const float mult_terms 	= grain_radius_term * /*Z_term */ n_H_term * T_term;

    /* abundance terms, depending on the key element for each grain species (normalised to solar ratio from Wiersma+09)*/
    /* !!! check consistency with yield factors, better to use an array of sola abundances defined in header? */


    if (p->id == 2149308301825) {
      message("mfrac :: | C %e | Si %e | O %e | Mg %e | Fe %e | Gra %e | Sil %e | Ide %e | Mgd %e | Fed %e",
    	      p->chemistry_data.metal_mass_fraction[chemistry_element_C],
    	      p->chemistry_data.metal_mass_fraction[chemistry_element_Si],
    	      p->chemistry_data.metal_mass_fraction[chemistry_element_O],
    	      p->chemistry_data.metal_mass_fraction[chemistry_element_Mg],
    	      p->chemistry_data.metal_mass_fraction[chemistry_element_Fe],
    	      p->chemistry_data.metal_mass_fraction[chemistry_element_Gra],
    	      p->chemistry_data.metal_mass_fraction[chemistry_element_Sil],
    	      p->chemistry_data.metal_mass_fraction[chemistry_element_Ide],
    	      p->chemistry_data.metal_mass_fraction[chemistry_element_Mgd],
    	      p->chemistry_data.metal_mass_fraction[chemistry_element_Fed]);
    }
    
    int di;
    float diffuse_fraction[5];

    /* Compute diffuse fraction of refractory elements for accretion */
    for (di = 0; di < 5; di++) {
      if (p->chemistry_data.metal_mass_fraction[dp->refractory_idx[di]] > 0){
	diffuse_fraction[di] = (p->chemistry_data.metal_mass_fraction[dp->refractory_idx[di]]/
				(p->chemistry_data.metal_mass_fraction[dp->dust_idx[di]]+
				 p->chemistry_data.metal_mass_fraction[dp->refractory_idx[di]]));
      }
      else {
	diffuse_fraction[di] =  0.;
      }
    }

    /* conversion factor converted to CGS from yr */
    const float tau_ref		= 6.31e15 * mult_terms;

    /* message("%f %f", fMgd, fFed);       */
    deltf_acc_Gra = dt_cgs * diffuse_fraction[0] / tau_ref; 
    deltf_acc_Ide = dt_cgs * diffuse_fraction[1] / tau_ref; 
    deltf_acc_Mgd = dt_cgs * diffuse_fraction[2] / tau_ref; 
    deltf_acc_Sil = dt_cgs * diffuse_fraction[3] / tau_ref; 
    deltf_acc_Fed = dt_cgs * diffuse_fraction[4] / tau_ref; 
    /* message("deltf_acc_Gra: %e", deltf_acc_Gra); */
  }
  else {
    deltf_acc_Gra = 0.; 
    deltf_acc_Sil = 0.; 
    deltf_acc_Ide = 0.; 
    deltf_acc_Mgd = 0.;
    deltf_acc_Fed = 0.; 
  }


  /* message("%e", max(1. + deltf_sput + deltf_SNII, 0.)); */
  /* compute total fractional change in mass for each grain species */
  double delta_mfrac_Gra = p->chemistry_data.metal_mass_fraction[chemistry_element_Gra] *
    (exp(deltf_sput + deltf_SNII + deltf_acc_Gra) - 1);
  double delta_mfrac_Sil = p->chemistry_data.metal_mass_fraction[chemistry_element_Sil] *
    (exp(deltf_sput + deltf_SNII + deltf_acc_Sil) - 1);
  double delta_mfrac_Ide = p->chemistry_data.metal_mass_fraction[chemistry_element_Ide] *
    (exp(deltf_sput + deltf_SNII + deltf_acc_Ide) - 1);
  double delta_mfrac_Mgd = p->chemistry_data.metal_mass_fraction[chemistry_element_Mgd] *
    (exp(deltf_sput + deltf_SNII + deltf_acc_Mgd) - 1);
  double delta_mfrac_Fed = p->chemistry_data.metal_mass_fraction[chemistry_element_Fed] *
    (exp(deltf_sput + deltf_SNII + deltf_acc_Fed) - 1);
  /* if (p->id == 3426890137225) { */
  /*   message("delta mfrac :: | Gra %e | Sil %e | Ide %e | Mgd %e | Fed %e | sput_tau %e | snii tau %e", */
  /* 	    delta_mfrac_Gra, delta_mfrac_Sil, delta_mfrac_Ide, delta_mfrac_Mgd, delta_mfrac_Fed, */
  /* 	    deltf_sput, deltf_SNII); */
  /* } */

  /* Message("%e %e %e", deltf_sput, deltf_SNII, deltf_acc_Gra); */
  /* message("%e %e", exp(deltf_sput+deltf_SNII+deltf_acc_Gra)-1, deltf_sput+deltf_SNII+deltf_acc_Gra); */
  /* if ((1-dp->with_accretion) && (delta_mfrac_Gra )){ */
  /*   message("%e %e", deltf_sput, deltf_SNII); */
  /* } */
  /* message("%e", (1-exp(deltf_sput + deltf_SNII + deltf_acc_Ide))); */
  /* ensure mass growth in dust does not exceed that of current mass in compositional element */

  /* delta_mfrac_Gra = min(delta_mfrac_Gra, */
  /* 			p->chemistry_data.metal_mass_fraction[chemistry_element_C]); */
  /* delta_mfrac_Ide = min(delta_mfrac_Ide, */
  /* 			p->chemistry_data.metal_mass_fraction[chemistry_element_Fe]); */
  /* delta_mfrac_Sil = min(delta_mfrac_Sil, */
  /* 			p->chemistry_data.metal_mass_fraction[chemistry_element_Mg]); */
  /* delta_mfrac_Mgd = min(delta_mfrac_Mgd, */
  /* 			p->chemistry_data.metal_mass_fraction[chemistry_element_Mg]); */
  /* delta_mfrac_Fed = min(delta_mfrac_Fed, */
  /* 			p->chemistry_data.metal_mass_fraction[chemistry_element_Fe]); */

  /* apply mass fraction changes to elements and grains */
  int i;
  for ( i = 0; i < chemistry_element_count; i++ ) {
    /* if (p->id == 3426890137225) { */
    /*   message("Element #%d mfrac: %e delta Element mfrac: %e", i, */
    /* 	      p->chemistry_data.metal_mass_fraction[i], */
    /* 	      (dp->comp_Gra[i]*delta_mfrac_Gra+ */
    /* 	       dp->comp_Sil[i]*delta_mfrac_Sil+ */
    /* 	       dp->comp_Ide[i]*delta_mfrac_Ide+ */
    /* 	       dp->comp_Mgd[i]*delta_mfrac_Mgd+ */
    /* 	       dp->comp_Fed[i]*delta_mfrac_Fed)); */
    /* } */

    p->chemistry_data.metal_mass_fraction[i] += dp->comp_Gra[i]*delta_mfrac_Gra;   
    p->chemistry_data.metal_mass_fraction[i] += dp->comp_Sil[i]*delta_mfrac_Sil;   
    p->chemistry_data.metal_mass_fraction[i] += dp->comp_Ide[i]*delta_mfrac_Ide;    
    p->chemistry_data.metal_mass_fraction[i] += dp->comp_Mgd[i]*delta_mfrac_Mgd;   
    p->chemistry_data.metal_mass_fraction[i] += dp->comp_Fed[i]*delta_mfrac_Fed;    

    /* if (p->id == 3426890137225) { */
    /*   message("Element #%d mfrac after: %e", i, */
    /* 	      p->chemistry_data.metal_mass_fraction[i]); */
    /* } */


    /* if (p->chemistry_data.metal_mass_fraction[i] < 0) { */
    if (0) {
      message("%d: %e ", i, 
	      p->chemistry_data.metal_mass_fraction[i]);
      message("Gra: %e, Sil %e, Ide %e, Mgd %e, Fed %e", 
	      delta_mfrac_Gra,delta_mfrac_Sil, delta_mfrac_Ide,
	      delta_mfrac_Mgd,delta_mfrac_Fed);
    /* message("factors: %e %e %e ", deltf_sput, deltf_SNII, */
    /* 	    deltf_acc_Gra+deltf_acc_Sil+deltf_acc_Ide); */   
      /* exit(0); */
    }

    
  }
  
  /* apply mass fraction changes to total metallicity */
  p->chemistry_data.metal_mass_fraction_total -= (delta_mfrac_Gra + delta_mfrac_Sil + delta_mfrac_Ide);
  
}


void dustevo_sputter_part(const struct phys_const *phys_const,
                       const struct unit_system *us,
                       const struct cosmology *cosmo,
                       const struct hydro_props *hydro_properties,
                       const struct entropy_floor_properties *floor_props,
                       const struct cooling_function_data *cooling,
		       const struct dustevo_props *dp,
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

  /* set subgrid particle properties if using */
  float rho, temp; 

  if (dp->with_subgrid_props) {
    rho = p->cooling_data.subgrid_dens; //xp->tracers_data.subgrid_dens;
    temp = p->cooling_data.subgrid_temp; //xp->tracers_data.subgrid_temp;
  }
  else {
    rho = hydro_get_physical_density(p, cosmo);
    temp = cooling_get_temperature(phys_const, hydro_properties, us, cosmo,
  					       cooling, p, xp);
  }

  /* Convert Hydrogen mass fraction into Hydrogen number density */
  const float XH = metal_fraction[chemistry_element_H];
  const double n_H = rho * XH / phys_const->const_proton_mass;
  const double n_H_cgs = n_H * units_cgs_conversion_factor(us, UNIT_CONV_NUMBER_DENSITY);
  const double dt_cgs = dt * units_cgs_conversion_factor(us, UNIT_CONV_TIME);

  /* surface area of a dust grain, assuming 0.1 micron radius */
  const double drdt_grain_cgs = -3.2e-18 * n_H_cgs * pow(pow(2e6/temp, 2.5) + 1, -1);
  const double dfdt_cgs = drdt_grain_cgs * dfdr_cgs;
  const double dfrac = dfdt_cgs * dt_cgs;

  /* p->chemistry_data.metal_mass_fraction[chemistry_element_Gra] *= max(1. + dfrac, 0.); */
  /* p->chemistry_data.metal_mass_fraction[chemistry_element_Sil] *= max(1. + dfrac, 0.); */
  /* p->chemistry_data.metal_mass_fraction[chemistry_element_Ide] *= max(1. + dfrac, 0.); */
  
  /* if (xp->sf_data.SFR > 0.) { */
  if (1<0){
  message("Sputter: %e", dfrac);
  }
}

void dustevo_accretion_part(const struct phys_const *phys_const,
                       const struct unit_system *us,
                       const struct cosmology *cosmo,
                       const struct hydro_props *hydro_properties,
                       const struct entropy_floor_properties *floor_props,
                       const struct cooling_function_data *cooling,
		       const struct dustevo_props *dp,
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

  /* set subgrid particle properties if using */
  float rho, temp; 

  if (dp->with_subgrid_props) {
    rho = p->cooling_data.subgrid_dens; //xp->tracers_data.subgrid_dens;
    temp = p->cooling_data.subgrid_temp; //xp->tracers_data.subgrid_temp;
  }
  else {
    rho = hydro_get_physical_density(p, cosmo);
    temp = cooling_get_temperature(phys_const, hydro_properties, us, cosmo,
  					       cooling, p, xp);
  }


  /* Convert Hydrogen mass fraction into Hydrogen number density */
  const double n_H = rho * X / phys_const->const_proton_mass;
  const double n_H_cgs = n_H * units_cgs_conversion_factor(us, UNIT_CONV_NUMBER_DENSITY);


  /* get other relevant particle properties */
  const float Z = p->chemistry_data.metal_mass_fraction_total;

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
  if (0 > 1) message("%e %e", temp, n_H_cgs); 
  if (0 > 1) message("%e %e", n_H_cgs, temp); 

  /* message("%e", 1./mult_terms);  */

  /* p->chemistry_data.metal_mass_fraction[chemistry_element_Gra] *= 1. + growfrac_Gra; */
  /* p->chemistry_data.metal_mass_fraction[chemistry_element_Sil] *= 1. + growfrac_Sil; */
  /* p->chemistry_data.metal_mass_fraction[chemistry_element_Ide] *= 1. + growfrac_Ide; */
  
  if (1<0){
    message("accretion: %e", growfrac_Gra+growfrac_Sil+growfrac_Ide);
  }
  /* exit(0); */
}

void dustevo_SNIIdestruction_part(const struct phys_const *phys_const,
				  const struct unit_system *us,
				  const struct cosmology *cosmo,
				  const struct hydro_props *hydro_properties,
				  const struct entropy_floor_properties *floor_props,
				  const struct cooling_function_data *cooling,
				  const struct dustevo_props *dp,
				  struct part *restrict p, struct xpart *restrict xp,
				  const float dt, const float dt_therm, const double time){
 
  const float X = p->chemistry_data.metal_mass_fraction[chemistry_element_H];
  const float fGra = p->chemistry_data.metal_mass_fraction[chemistry_element_Gra];
  const float fSil = p->chemistry_data.metal_mass_fraction[chemistry_element_Sil];
  const float fIde = p->chemistry_data.metal_mass_fraction[chemistry_element_Ide];

  /* no dust, no SNII destruction */
  if (fGra+fSil+fIde == 0.) return;

  /* no SF, no SNII */
  if (xp->sf_data.SFR <= 0) return;

  /* constants */ 
  const float eps_eff = 0.1;
  const double unit_mass_cgs = units_cgs_conversion_factor(us, UNIT_CONV_MASS);
  const double unit_time_cgs = units_cgs_conversion_factor(us, UNIT_CONV_TIME);
  const double dt_cgs = dt * unit_time_cgs;


  /* set subgrid particle properties if using */
  float rho, temp; 
  if (dp->with_subgrid_props) {
    rho = p->cooling_data.subgrid_dens; //xp->tracers_data.subgrid_dens;
    temp = p->cooling_data.subgrid_temp; //xp->tracers_data.subgrid_temp;
  }
  else {
    rho = hydro_get_physical_density(p, cosmo);
    temp = cooling_get_temperature(phys_const, hydro_properties, us, cosmo,
  					       cooling, p, xp);
  }

  /* Convert Hydrogen mass fraction into Hydrogen number density */
  const double n_H = rho * X / phys_const->const_proton_mass;
  const double n_H_cgs = n_H * units_cgs_conversion_factor(us, UNIT_CONV_NUMBER_DENSITY);
  /* const double n_H_sg_cgs = n_H_cgs * (rho_sg/rho); */

  /* get other relevant particle properties */
  const float Z = p->chemistry_data.metal_mass_fraction_total;

  /* SNII rate assuming constant SFR, in cgs units */
  const float rate_SNII = xp->sf_data.SFR * dp->specific_numSNII / unit_time_cgs;
  
  /* Yamasawa et. al (2011): the cgs mass of material swept up in SNII shocks */
  const float m_swept = 3.052e+36 * pow(n_H_cgs, -0.202)* pow((Z/0.0127) + 0.039, -0.289);

  /* destruction timescale of dust in cgs */  
  const float tau_dest = p->mass * unit_mass_cgs / (eps_eff*m_swept * rate_SNII);

  /* fractional change in dust due to dust destruction */
  const float dfrac = dt_cgs/tau_dest;

  /* p->chemistry_data.metal_mass_fraction[chemistry_element_Gra] *= max(1. - dfrac, 0); */
  /* p->chemistry_data.metal_mass_fraction[chemistry_element_Sil] *= max(1. - dfrac, 0); */
  /* p->chemistry_data.metal_mass_fraction[chemistry_element_Ide] *= max(1. - dfrac, 0); */

  /* message("SNII : %e", dfrac); */
  if (1<0){ /* (xp->sf_data.SFR > 0.) { */
    message("SNII: %e", m_swept);
    message("SNII: %e", xp->sf_data.SFR);
    message("SNII: %e", dfrac);
    message("SNII: %e", temp);
    exit(0);
  }

  /* message("Exiting..."); */
  /* exit(0); */
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

void redistribute_dust_to_element_abundances(struct spart* sp, 
					     const struct dustevo_props *dp){

  /* return mass from dust to individual elements according to assumed grain composition */
  if (dp->model_type == 0){
    /* message("redistributing according to dust model %d...", dp->model_type); */
    sp->chemistry_data.metal_mass_fraction[chemistry_element_C] += 
      sp->chemistry_data.metal_mass_fraction[chemistry_element_Gra];

    sp->chemistry_data.metal_mass_fraction[chemistry_element_Fe] += 
      sp->chemistry_data.metal_mass_fraction[chemistry_element_Ide];

    sp->chemistry_data.metal_mass_fraction[chemistry_element_O] += 
      0.3715755*sp->chemistry_data.metal_mass_fraction[chemistry_element_Sil];
    sp->chemistry_data.metal_mass_fraction[chemistry_element_Mg] += 
      0.1411169*sp->chemistry_data.metal_mass_fraction[chemistry_element_Sil]; 
    sp->chemistry_data.metal_mass_fraction[chemistry_element_Si] += 
      0.1630668*sp->chemistry_data.metal_mass_fraction[chemistry_element_Sil];
    sp->chemistry_data.metal_mass_fraction[chemistry_element_Fe] += 
      0.3242408*sp->chemistry_data.metal_mass_fraction[chemistry_element_Sil];
  }

  if (dp->model_type == 1){
    /* message("redistributing according to dust model %d...", dp->model_type); */
    sp->chemistry_data.metal_mass_fraction[chemistry_element_C] += 
      sp->chemistry_data.metal_mass_fraction[chemistry_element_Gra];

    sp->chemistry_data.metal_mass_fraction[chemistry_element_Si] += 
      sp->chemistry_data.metal_mass_fraction[chemistry_element_Sil];

    sp->chemistry_data.metal_mass_fraction[chemistry_element_O] += 
      sp->chemistry_data.metal_mass_fraction[chemistry_element_Ide];

    sp->chemistry_data.metal_mass_fraction[chemistry_element_Mg] += 
      sp->chemistry_data.metal_mass_fraction[chemistry_element_Mgd];

    sp->chemistry_data.metal_mass_fraction[chemistry_element_Fe] += 
      sp->chemistry_data.metal_mass_fraction[chemistry_element_Fed];

    /* account for returned nebular metals from extra dust channels in total metal mass fraction */
    sp->chemistry_data.metal_mass_fraction_total +=
      sp->chemistry_data.metal_mass_fraction[chemistry_element_Mgd] +
      sp->chemistry_data.metal_mass_fraction[chemistry_element_Fed];
    
    /* set mass fraction of extra dust channels inside stellar atmospheres to zero */
    sp->chemistry_data.metal_mass_fraction[chemistry_element_Mgd] = 0.;
    sp->chemistry_data.metal_mass_fraction[chemistry_element_Fed] = 0.;
    
  }

  /* account for returned nebular metals from universal dust types in total metal mass fraction */
  sp->chemistry_data.metal_mass_fraction_total +=
    sp->chemistry_data.metal_mass_fraction[chemistry_element_Gra] +
    sp->chemistry_data.metal_mass_fraction[chemistry_element_Sil] +
    sp->chemistry_data.metal_mass_fraction[chemistry_element_Ide];

  /* set mass fraction of universal dust grains inside stellar atmospheres to zero */
  sp->chemistry_data.metal_mass_fraction[chemistry_element_Gra] = 0.;
  sp->chemistry_data.metal_mass_fraction[chemistry_element_Sil] = 0.;
  sp->chemistry_data.metal_mass_fraction[chemistry_element_Ide] = 0.;
}
