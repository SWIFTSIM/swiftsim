/* Local includes. */
#include "chemistry.h"
#include "cooling.h"
#include "dust_properties.h"
#include "dust_yield_tables.h"
#include "dust.h"
#include <string.h>


/**
 * @brief redistribute any dust mass back to element abundances
 * on star particle formation according to dust composition, to 
 * represent astration 
 *
 * @param p The gas particles.  
 * @param dp Global dust parameters for initialisation.
 */
void redistribute_dust_masses(struct spart* sp,
			      const struct part* p, 
			      const struct dustevo_props *dp) {
  
  /** 
   * iterate through grain species and element types and
   * redistribute dust abundances to element abundances,
   * according to composition array
   */

  float grain_mass;
  int eldx;
  for (int grain = 0; grain < grain_species_count; grain++) {
    grain_mass =  p->dust_data.grain_mass_fraction[grain];
    for (int elem = 0; elem < dp->grain_element_count[grain]; elem++) {
      eldx = dp->grain_element_indices[grain][elem];
      /* put fraction of grain mass back into constituent element (astration) */
      sp->chemistry_data.metal_mass_fraction[eldx] +=
	grain_mass * dp->grain_element_mfrac[grain][elem]; 
    }
  }    
}

/**
 * @brief Prints the dust evolution model to stdout.
 *
 * @param dust #dustevo_props struct.
 */
void dustevo_print_backend(const struct dustevo_props *dp) {
    message("Running with a Trayford et al. (2020) dust evolution model.");
}

 /**
 * @brief initialise structure housing global dust parametrisation.
 * In particular, flags and values set in the parameter file, 
 * and any hard-coded properties
 * @brief Prints the dust evolution model to stdout.
 *
 * @param dp Global dust parameters for initialisation
 * @param params The parsed parameter file.
 * @param phys_const The physical constants in internal units.
 * @param us The current internal system of units.
 * @param dust #dustevo_props struct.
 */
void dustevo_props_init_backend(struct dustevo_props* dp,
				struct swift_params* params,
				struct feedback_props* fp,
				struct cooling_function_data* cooling,
				const struct phys_const* phys_const,
				const struct unit_system* us) {

  /**
   * initialise structure housing global dust parametrisation.
   * First, set global parameters: 
   *
   * initialise dustevo_props struct object and store:
   * - Parameter file flags (e.g. with_accretion, with_sputtering)
   * - Parameter file values (e.g diffusion boost, clumping factor)
   * - Hard coded values (e.g. grain radii, density, composition)
   *
   *
   * Then, using the yield tables from feedback_props, compute 
   * and set dust yields for each of the  dust species 
   *
   * Finally apply corections to the cooling and heating rates
   * to remove implicit dust, via the scale_out_table_depletion
   * function.
   **/  
  
  /* read some parameters */
  /* operation modes */
  dp->pair_to_cooling = 
    parser_get_opt_param_int(params, "DustEvolution:pair_to_cooling", 0);
  dp->with_subgrid_props = 
    parser_get_opt_param_int(params, "DustEvolution:use_subgrid_props", 1);
  dp->with_sputtering = 
    parser_get_opt_param_int(params, "DustEvolution:use_sputtering", 1);
  dp->with_SNII_destruction = 
    parser_get_opt_param_int(params, "DustEvolution:use_SNII_destruction", 1);
  dp->with_accretion = 
    parser_get_opt_param_int(params, "DustEvolution:use_accretion", 1);

  /* initial abundances */
  dp->initial_grain_mass_fraction[grain_species_graphite] = 
    parser_get_opt_param_float(params, "DustEvolution:initial_abundance_graphite", 0.);
  dp->initial_grain_mass_fraction[grain_species_silicate] =
    parser_get_opt_param_float(params, "DustEvolution:initial_abundance_silicate", 0.);

  /* global parameters */
  dp->clumping_factor = 
    parser_get_opt_param_float(params, "DustEvolution:clumping_factor", 10.);
  dp->diffusion_rate_boost = 
    parser_get_opt_param_float(params, "DustEvolution:diffusion_boost_factor",
			       1.0);
  dp->nu = 
    parser_get_opt_param_float(params, "DustEvolution:silicate_fe_grain_fraction",
			       0.5);

  parser_get_opt_param_string(params, "DustEvolution:dust_yields_path", 
			      dp->AGB_dyield_path, "./dust_yields");
  
  /* set assumed abundance patterns to Wiersma et al (2009a) */
  memset(dp->abundance_pattern, 0, sizeof dp->abundance_pattern);
  dp->abundance_pattern[chemistry_element_H] = 7.0649785e-01;
  dp->abundance_pattern[chemistry_element_He] = 2.8055534e-01;
  dp->abundance_pattern[chemistry_element_C] = 2.0665436e-03;
  dp->abundance_pattern[chemistry_element_N] = 8.3562563e-04;
  dp->abundance_pattern[chemistry_element_O] = 5.4926244e-03;
  dp->abundance_pattern[chemistry_element_Ne] = 1.4144605e-03;
  dp->abundance_pattern[chemistry_element_Mg] = 5.9070642e-04;
  dp->abundance_pattern[chemistry_element_Si] = 6.8258739e-04;
  dp->abundance_pattern[chemistry_element_Fe] = 1.1032152e-03;

  dp->solar_metallicity = 0.0127;

  /* set element atomic weights */
  memset(dp->atomic_weight, 0, sizeof dp->atomic_weight);
  dp->atomic_weight[chemistry_element_H] = 1.0079;
  dp->atomic_weight[chemistry_element_He] = 4.0026;
  dp->atomic_weight[chemistry_element_C] = 12.0107;
  dp->atomic_weight[chemistry_element_N] = 14.0067;
  dp->atomic_weight[chemistry_element_O] = 15.9994;
  dp->atomic_weight[chemistry_element_Ne] = 20.1797;
  dp->atomic_weight[chemistry_element_Mg] = 24.305;
  dp->atomic_weight[chemistry_element_Si] = 28.0855;
  dp->atomic_weight[chemistry_element_Fe] = 55.845;

  /* set element count contributing to each grain */
  dp->grain_element_count[grain_species_graphite] = 1;
  dp->grain_element_count[grain_species_silicate] = 4;

  /* effective grain molecular weights in amu */
  float Agra = dp->atomic_weight[chemistry_element_C];
  float Asil = 2* dp->nu * dp->atomic_weight[chemistry_element_Fe] +
    (2 - 2*dp->nu) * dp->atomic_weight[chemistry_element_Mg] +
    dp->atomic_weight[chemistry_element_Si] +
    4*dp->atomic_weight[chemistry_element_O];


  /* set condensation fractions */
  memset(dp->condensation_frac, 0, sizeof dp->condensation_frac);
  dp->condensation_frac[grain_species_graphite] = 0.15; 
  dp->condensation_frac[grain_species_silicate] = 3.5e-4;

  /* set grain composition */
  float gcomp1[1] = {dp->atomic_weight[chemistry_element_C]/Agra};
  float gcomp2[4] = {4*dp->atomic_weight[chemistry_element_O]/Asil,
		     (2-2*dp->nu)*dp->atomic_weight[chemistry_element_Mg]/Asil,
  		     dp->atomic_weight[chemistry_element_Si]/Asil,
  		     2*dp->nu*dp->atomic_weight[chemistry_element_Fe]/Asil};

  dp->grain_element_mfrac[grain_species_graphite] = gcomp1;
  dp->grain_element_mfrac[grain_species_silicate] = gcomp2;

  /* allocate memory */
  dp->grain_element_mfrac[grain_species_graphite] = 
  malloc(dp->grain_element_count[grain_species_graphite] * sizeof(float)); 
  dp->grain_element_mfrac[grain_species_silicate] = 
  malloc(dp->grain_element_count[grain_species_silicate] * sizeof(float));  

  /* deep copy grain arrays */
  memcpy(dp->grain_element_mfrac[grain_species_graphite],
	 &gcomp1, dp->grain_element_count[grain_species_graphite] * sizeof(float));
  memcpy(dp->grain_element_mfrac[grain_species_silicate],
	 &gcomp2, dp->grain_element_count[grain_species_silicate] * sizeof(float));

  /* set chemistry indices composition */
  int gidx1[1] = {chemistry_element_C};
  int gidx2[4] = {chemistry_element_O,
  		  chemistry_element_Mg,
  		  chemistry_element_Si,
  		  chemistry_element_Fe};


  /* allocate memory */
  dp->grain_element_indices[grain_species_graphite] = 
  malloc(dp->grain_element_count[grain_species_graphite] * sizeof(int)); 
  dp->grain_element_indices[grain_species_silicate] = 
  malloc(dp->grain_element_count[grain_species_silicate] * sizeof(int));  

  /* deep copy grain arrays */
  memcpy(dp->grain_element_indices[grain_species_graphite],
	 &gidx1, dp->grain_element_count[grain_species_graphite] * sizeof(int));
  memcpy(dp->grain_element_indices[grain_species_silicate],
	 &gidx2, dp->grain_element_count[grain_species_silicate] * sizeof(int));

  /* set element count contributing to each grain */
  dp->accretion_coeff[grain_species_graphite] = 3.132e15;
  dp->accretion_coeff[grain_species_silicate] = 5.677e15;

  /* number of SNII per unit stellar mass */
  dp->specific_numSNII = 7.039463e-3 / phys_const->const_solar_mass;

  /** NOTE 1: total metallicity yields untouched here, so Z represents the conserved dust + gas phase metals **/

  /** NOTE 2: only the IMF resampled tables are modified in fp, while plain .yield arrays are unchanged (the 
   *  original yield tables are only used to compute these and are already modified via the SNII yield factors) **/

  initialise_dyield_tables(fp, dp);
  compute_AGB_dyield(fp, dp);
  compute_SNII_dyield(fp, dp);
  //print_dyield_tables(fp, dp);
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
  int eldx;
  
  /* metal mass fraction in dust before evolution step */
  float D_pre = 0.;
  for (int grain = 0; grain < grain_species_count; grain++) {
    D_pre += p->dust_data.grain_mass_fraction[grain];
  }
  
  /* no dust, no accretion or destruction. return. */
  if (D_pre <= 0.) return;

  /* set temperature and density corresponding to either subgrid or hydro */
  float rho, temp;

  if (dp->with_subgrid_props) {
    rho = p->cooling_data.subgrid_dens;
    temp = p->cooling_data.subgrid_temp;
  }
  else {
    rho = hydro_get_physical_density(p, cosmo);
    temp = cooling_get_temperature(phys_const, hydro_properties, us, cosmo,
  					       cooling, p, xp);
  }

  /* --------------- Common Parameters --------------- */

  /* grain radius in cgs, hard coded for now, 0.1 micron. */
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
    const float m_swept = 3.052e+36 * pow(n_H_cgs, -0.202) * 
      pow((Z/dp->solar_metallicity) + 0.039, -0.298);

    /* destruction timescale of dust in cgs */
    const float tau_dest = p->mass * unit_mass_cgs / (eps_eff*m_swept * rate_SNII);

    /* fractional change in dust due to dust destruction */
    deltf_SNII = -dt_cgs/tau_dest;
  }
  else {
    deltf_SNII = 0.;
  }

  /* ------------------------- Accretion ------------------------- */
  
  double deltf_acc[grain_species_count] = {0.};

  if (dp->with_accretion){
    /* accretion timescale multiplicative terms  (Hirashita & Voshchinnikov 2013),
       implicitly assuming S_acc = 0.3 */
    const float grain_radius_term = grain_radius_cgs / 1e-5;
    const float n_H_term	  = 1.e3 / n_H_cgs;
    const float T_term		  = pow(10./temp, 0.5);
    const float mult_terms 	  = grain_radius_term * n_H_term * T_term;
    float min_abundance_ratio = FLT_MAX;
    float tau_acc;

    for (int grain = 0; grain < grain_species_count; grain++) {
      for (int elem = 0; elem < dp->grain_element_count[grain]; elem++) {
	eldx = dp->grain_element_indices[grain][elem];
	/* find the minimum solar abundance ratio among constituent elements 
	   indicating element that bottlenecks accretion of effective grain */
	min_abundance_ratio = min(min_abundance_ratio, 
	  p->chemistry_data.metal_mass_fraction[eldx] / dp->abundance_pattern[eldx]);
      }
      tau_acc = dp->accretion_coeff[grain] * mult_terms / min_abundance_ratio;
      deltf_acc[grain] = dt_cgs  * dp->clumping_factor / (tau_acc);
    }
  }
  
  /* ------------------------- Combine ------------------------- */
  
  float delta_dfrac;
  float D_post = 0.;

  for (int grain = 0; grain < grain_species_count; grain++) {
    delta_dfrac = p->dust_data.grain_mass_fraction[grain] *
      (exp(deltf_sput + deltf_SNII + deltf_acc[grain]) - 1.);
    /* limit grain growth */
    for (int elem = 0; elem < dp->grain_element_count[grain]; elem++) {
      eldx = dp->grain_element_indices[grain][elem];
      /* grain growth limited by available mass in a consitituent element */
      delta_dfrac = min(delta_dfrac,
			p->chemistry_data.metal_mass_fraction[eldx] /
			dp->grain_element_mfrac[grain][elem]);
      /* message("eldx %d mfrac %e", eldx, dp->grain_element_mfrac[grain][elem]); */
    }

    /* update particle grain mass fractions */
    p->dust_data.grain_mass_fraction[grain] += delta_dfrac;
 
    for (int elem = 0; elem < dp->grain_element_count[grain]; elem++) {
      eldx = dp->grain_element_indices[grain][elem];
      /* update constituent element abundances */
      p->chemistry_data.metal_mass_fraction[eldx] -= delta_dfrac *
	dp->grain_element_mfrac[grain][elem];
    }

    D_post += p->dust_data.grain_mass_fraction[grain];
  }
  if (0) {
    message("delta dust: %e, T %e rho %e", (D_post-D_pre)/D_pre, temp, rho);
  }
}
