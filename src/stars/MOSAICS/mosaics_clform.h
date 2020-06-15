/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2019 Joel Pfeffer (j.l.pfeffer@ljmu.ac.uk)
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
#ifndef SWIFT_MOSAICS_CLFORM_H
#define SWIFT_MOSAICS_CLFORM_H

#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>

/* Local includes */
#include "star_formation.h"
#include "physical_constants.h"
#include "dsyevj3.h"
#include "mosaics_cfelocal.h"

/**
 * @brief Generate a uniform random number
 */
__attribute__((always_inline)) INLINE static double get_rand(
    gsl_rng* random_generator) {

  return gsl_rng_uniform(random_generator);
}

/**
 * @brief Normalisation of Schechter function for probability function
 */
__attribute__((always_inline)) INLINE static double schechternorm(
    double mmin, double mmax, double slope, double mstar) {

  /* integration resolution */
  const int nintsteps = 10000;

  /* starting mass */
  double mc = mmin;
  double norm = 0.f;

  /* Constant integration step in log space */
  const double dm = pow(mmax/mmin, 1.f/(double)nintsteps) - 1.f;

  /* integrate to determine norm (total number) */
  for (int n=0; n<nintsteps; n++) {
    /* dm */
    double dmc = mc*dm;

    /* Leapfrog */
    norm += pow(mc,-slope) * exp(-mc/mstar) * dmc -
        0.5 * (pow(mc,-slope) / mstar + slope * pow(mc,-slope-1.f)) *
        exp(-mc/mstar) * (dmc*dmc);

    mc += dmc;
  }

  return norm;
}

/**
 * @brief Mean mass of Schechter function
 *
 * Expected cluster mass \int_{mmin}^{mmax} M N(M) dM
 */
__attribute__((always_inline)) INLINE static double meanmass(
    double mmin, double mmax, double slope, double mstar, double norm) {

  /* integration resolution*/
  const int nintsteps = 10000;

  /* starting mass */
  double mc = mmin;
  double meanm = 0.f;

  /* Constant integration step in log space */
  const double dm = pow(mmax/mmin, 1.f/(double)nintsteps) - 1.f;

  /* integrate to determine norm (total number) */
  for (int n=0; n<nintsteps; n++) {
    /* dm */
    double dmc = mc*dm;

    /* Leapfrog */
    meanm += pow(mc,-slope+1.f) * exp(-mc/mstar) * dmc -
        0.5 * (pow(mc,-slope+1.f) / mstar + (slope-1.f) * pow(mc,-slope)) *
        exp(-mc/mstar) * (dmc*dmc);

    mc += dmc;
  }

  return meanm/norm;
}

/**
 * @brief Draw masses from Schechter function
 */
__attribute__((always_inline)) INLINE static double schechtergen (
    double mmin, double mmax, double slope, double mstar, double norm,
    gsl_rng* random_generator) {

  /* integration resolution */
  const int nintsteps = 10000;

  const double x = get_rand(random_generator);

  /* starting mass */
  double mc = mmin;
  double mcnorm = 0.f;

  /* Constant integration step in log space */
  const double dm = pow(mmax/mmin, 1.f/(double)nintsteps)-1.f;

  /* integrate until we reach the random number */
  double diff = 0.f;
  double dmc = 0.f;
  while (mcnorm/norm < x) {
    /* dm */
    dmc = mc*dm;

    /* leapfrog */
    diff = pow(mc,-slope) * exp(-mc/mstar) * dmc - 
        0.5 * (pow(mc,-slope) / mstar + slope * pow(mc,-slope-1.f)) * 
        exp(-mc/mstar) * (dmc*dmc);

    mcnorm += diff;
    mc += dmc;
  }

  /* Interpolate between integration steps to correct mass */
  /* Previous step */
  const double x1 = (mcnorm - diff) / norm;
  /* Present step */
  const double x2 = mcnorm / norm;
  const double frac = (x-x1) / (x2-x1);

  /* interpolated mass */
  return pow(mc-dmc, 1.f-frac) * pow(mc, frac);
}

/**
 * @brief Draw masses from power-law 
 */
__attribute__((always_inline)) INLINE static double powerlawgen(
    double mmin, double mmax, double slope, gsl_rng* random_generator) {

  double x = get_rand(random_generator);

  return pow(x * (pow(mmax,1.f-slope) - pow(mmin,1.f-slope)) + 
      pow(mmin,1.f-slope), 1.f/(1.f-slope));
}

/**
 * @brief Returns the value ln[Gamma( xx )] for xx > 0.
 *
 * From Numerical Recipes in C 
 */
__attribute__((always_inline)) INLINE static double gammln(double xx) {

    /* Internal arithmetic will be done in double precision, a nicety that
     * you can omit if five-figure accuracy is good enough. */
    double x,y,tmp,ser;
    static double cof[6]={76.18009172947146,-86.50532032941677,
        24.01409824083091,-1.231739572450155,
        0.1208650973866179e-2,-0.5395239384953e-5};
    int j;

    y=x=xx;
    tmp=x+5.5;
    tmp -= (x+0.5)*log(tmp);
    ser=1.000000000190015;
    for (j=0;j<=5;j++) ser += cof[j]/++y;
    return -tmp+log(2.5066282746310005*ser/x);
}

/**
 * @brief Draw a random value from Poisson distribution
 *
 * Returns as a floating-point number an integer value that is a random deviate 
 * drawn from a Poisson distribution of mean xm.
 * From Numerical Recipes in C
 */
__attribute__((always_inline)) INLINE static double poidev(
    double xm, gsl_rng* random_generator) {

    static double sq,alxm,g,oldm=(-1.0); /* oldm is a flag for whether xm has changed      since last call. */
    double em,t,y;
    if (xm < 12.0) {
        /* Use direct method. */
        if (xm != oldm) {
            oldm=xm;
            g=exp(-xm); /* If xm is new, compute the exponential. */
        }
        em = -1;
        t=1.0;
        do {
            /* Instead of adding exponential deviates it is equivalent to multiply 
               uniform deviates. We never actually have to take the log, merely 
               compare to the pre-computed exponential. 
             */
            ++em;
            t *= get_rand(random_generator);
        } while (t > g);
    } else {
        /* Use rejection method. */
        if (xm != oldm) {
            /* If xm has changed since the last call, then precompute some functions that  occur below. */
            oldm=xm;
            sq=sqrt(2.0*xm);
            alxm=log(xm);
            g=xm*alxm-gammln(xm+1.0);
            /* The function gammln is the natural log of the gamma function, as given in 6.1. */
        }
        do {
            do {
                /* y is a deviate from a Lorentzian comparison function. */
                y=tan(M_PI*get_rand(random_generator));
                em=sq*y+xm; /* em is y, shifted and scaled. */
            } while (em < 0.0); /* Reject if in regime of zero probability. */
            em=floor(em); /* The trick for integer-valued distributions. */
            t=0.9*(1.0+y*y)*exp(em*alxm-gammln(em+1.0)-g);
            /* The ratio of the desired distribution to the comparison function; we accept or
            reject by comparing it to another uniform deviate. The factor 0.9 is chosen so
            that t never exceeds 1.
             */
        } while (get_rand(random_generator) > t);
    }
    return em;
}

/**
 * @brief Do the cluster formation for the mosaics subgrid star cluster model
 *
 * @param sp The particle to act upon
 * @param stars_properties Properties of the stars model.
 * @param sf_props the star formation law properties to use
 * @param phys_const the physical constants in internal units.
 */
__attribute__((always_inline)) INLINE static void mosaics_clform(
    struct spart* restrict sp, const struct engine* e,
    const int with_cosmology) {

#if !defined(STAR_FORMATION_NONE)

  const struct stars_props* props = e->stars_properties;
  const struct star_formation *sf_props = e->star_formation;
  const struct phys_const *phys_const = e->physical_constants;
  const struct cosmology *cosmo = e->cosmology;
  const struct unit_system *us = e->internal_units;

  /* Some constants */
  const double const_G = phys_const->const_newton_G;
  const double density_to_kgm3 =
      units_cgs_conversion_factor(us, UNIT_CONV_DENSITY) * 1e3;
  const double velocity_to_ms =
      units_cgs_conversion_factor(us, UNIT_CONV_VELOCITY) * 0.01;
  const double time_to_cgs = units_cgs_conversion_factor(us, UNIT_CONV_TIME);

  /* TODO for now just use the particle ID as seed */
  /* Set up the random number generator */
  gsl_rng *random_generator = gsl_rng_alloc(gsl_rng_ranlxd1);
  gsl_rng_set(random_generator, sp->id);

  /* Finish the velocity dispersion/gas fraction calculation */
  const float h_inv = 1.f / sp->hbirth;
  const float h_inv2 = h_inv * h_inv;

  if (sp->scount > 0) {
    sp->stars_rho *= h_inv * h_inv2;
    sp->stars_sigma_v2 *= (1.f/sp->stars_rho) * h_inv * h_inv2 * cosmo->a2_inv;
  }

  sp->fgas = sp->gas_mass_unweighted / 
      (sp->gas_mass_unweighted + sp->stars_mass_unweighted);

  const double star_vel_disp = sqrt(sp->stars_sigma_v2);
  const double gas_vel_disp = sqrt(sp->sf_data.birth_velocity_dispersion);

  /* -------- Get CFE -------- */
  /* In units of kg, m, s */

  /* The local gas density */
  double rholoc = sp->sf_data.birth_density * density_to_kgm3;

  /* The local gas velocity dispersion */
  double sigmaloc;
  if (props->subgrid_gas_vel_disp) {
    /* "Sub-particle turbulent velocity dispersion" */
    sigmaloc = sqrt(sp->birth_pressure / sp->sf_data.birth_density);
  } else {
    /* The resolved velocity dispersion */
    /* 1/sqrt(3) converts to 1D assuming isotropy */
    sigmaloc = gas_vel_disp / sqrt(3.f);
  }
  sigmaloc *= velocity_to_ms;

  /* Sound speed of cold ISM */
  double csloc;
  if (props->Fixedcs > 0) {
    csloc = props->Fixedcs;
  } else {
    csloc = sp->sound_speed_subgrid * velocity_to_ms;
  }

  /* Use a fixed value for the CFE? */
  if (props->FixedCFE > 0) {
    sp->CFE = props->FixedCFE;
  } else {
    /* Calculate CFE based on local conditions (Kruijssen 2012). units kg, m, s*/
    sp->CFE = f_cfelocal(rholoc, sigmaloc, csloc, props);
  }

  /* Get feedback timescale while we're here and convert to internal units */
  double tfb = feedback_timescale(rholoc, sigmaloc, csloc, props);
  tfb /= time_to_cgs;

  /* -------- Get Mcstar -------- */

  /*TODO this is all unneccessary if we have a power law ICMF */

  /* Gas surface density (Krumholz & McKee 2005) */
  double phi_P = 1.f;
  /*TODO work out the minimum number we need */
  /* Make sure we have enough stars for a resolved dispersion */
  if ( sp->scount > 5 ) {
    phi_P = 1.f + gas_vel_disp / star_vel_disp * (1.f / sp->fgas - 1.f);
  }

  const double total_pressure = 
      sp->sf_data.birth_density * sp->sf_data.birth_velocity_dispersion / 3.f;

  const double SigmaG = sqrt(2. * total_pressure / (M_PI * const_G * phi_P));

  /* Toomre mass via tidal tensors (Pfeffer+18) */

  /* Orbital frequencies */
  double kappa2;

  /* Are we using the smoothed version? And is there any stellar neighbours? */
  if (props->smoothed_orbit_frequencies && sp->scount > 0) {
    /* Finish the neighbour calculations */
    sp->Omega_birth *= (1.f/sp->stars_rho) * h_inv * h_inv2;
    sp->kappa_birth *= (1.f/sp->stars_rho) * h_inv * h_inv2;

  } else {
    /* No neighbours, so revert to value from own tensor */

    /* Temporary array for eigvec/val calculation
     * Note the lower elements are never referenced by dsyevj3 */
    double tensor[3][3];
    tensor[0][0] = sp->tidal_tensor[2][0];
    tensor[0][1] = sp->tidal_tensor[2][1];
    tensor[0][2] = sp->tidal_tensor[2][2];
    tensor[1][1] = sp->tidal_tensor[2][3];
    tensor[1][2] = sp->tidal_tensor[2][4];
    tensor[2][2] = sp->tidal_tensor[2][5];

    /* Get the eigenvalues */
    double tideval[3];
    dsyevj3(tensor, tideval);

    /* sort eigenvalues, tideval[2] (=lambda_1) is largest eigenvalue */
    sort3(tideval);

    /* Circular and epicyclic frequencies */
    double Omega2;
    if (props->Omega_is_lambda2) {
      /* Correct version (in principle) */
      Omega2 = fabs(-tideval[1]);
    } else {
      /* Version in E-MOSAICS */
      Omega2 = fabs(-tideval[0] - tideval[1] - tideval[2]) / 3.f;
    }
    kappa2 = fabs(3.f * Omega2 - tideval[2]);

    sp->Omega_birth = sqrt(Omega2);
    sp->kappa_birth = sqrt(kappa2);
  }

  //TODO Birth Omega and kappa should not have scale factors

  /* Factor out cosmology */
  kappa2 = sp->kappa_birth * sp->kappa_birth * cosmo->a3_inv;

  /* 4 * pi^5 * G^2 */
  const double MTconst =
      4.f * M_PI * M_PI * M_PI * M_PI * M_PI * const_G * const_G;
  sp->Toomre_mass = MTconst * SigmaG * SigmaG * SigmaG / (kappa2 * kappa2);

  /* Total collapse time of Toomre volume */
  const double tcollapse = sqrt(2. * M_PI / kappa2);

  /* Toomre collapse fraction (Reina-Campos & Kruijssen 2017) */
  sp->frac_collapse = fmin(1., tfb / tcollapse);
  sp->frac_collapse *= sp->frac_collapse; /* f^2 */
  sp->frac_collapse *= sp->frac_collapse; /* f^4 */

  /* The 'GMC' mass */
  const double M_collapse = sp->Toomre_mass * sp->frac_collapse;

  /* Exponential truncation of cluster mass function */
  if (props->SFE > 0) {

    sp->Mcstar = props->SFE * sp->CFE * M_collapse;

  } else {

    /* TODO untested model... */

#if defined(STAR_FORMATION_COLIBRE)
    double sfe_ff;
    switch (sf_props->SF_law) {
      case colibre_star_formation_schmidt_law:
        sfe_ff = sf_props->schmidt_law.sfe;
        break;
      default:
        sfe_ff = 0.01;
        break;
    }
#elif defined(STAR_FORMATION_GEAR)
    const double sfe_ff = sf_props->star_formation_efficiency;
#else
    /* No SFE defined for star formation model */
    /* Elmegreen 2002 */
    const double sfe_ff = 0.012;
#endif

    /* Free-fall time */
    const double tff = 
        sqrt(3.0 * M_PI / (32.0 * const_G * sp->sf_data.birth_subgrid_density));

    /* Integrated SFE. End of collapse defined by shortest of Toomre collapse
     * timescale and feedback timescale */
    const double SFE_int = sfe_ff * fmin(tcollapse, tfb) / tff;

    if ((SFE_int < 0.0) || (SFE_int > 1.0)) {
      error("Integrated SFE unphysical! SFE_int=%g, sfe_ff=%g, tff=%g,"
            "tfb=%g, tcoll=%g.",SFE_int, sfe_ff, tff, tfb, tcollapse);
    }

    sp->Mcstar = SFE_int * sp->CFE * M_collapse;
  }


  /* ----- Set up the star cluster population for this particle ----- */

  /* TODO Where's the best place to initialise everything? */
  sp->initial_num_clusters = 0;
  sp->initial_num_clusters_evo = 0;
  sp->num_clusters = 0;
  sp->initial_cluster_mass_total = 0.f;
  sp->initial_cluster_mass_evo = 0.f;
  sp->cluster_mass_total = 0.f;
  for (int i = 0; i < MOSAICS_MAX_CLUSTERS; i++) {
    //sp->clusters.id[i] = -1;
    sp->clusters.mass[i] = 0.f;
    sp->clusters.initial_mass[i] = 0.f;
    sp->clusters.dmevap[i] = 0.f;
    sp->clusters.dmshock[i] = 0.f;
    sp->clusters.disruption_time[i] = -1.f;
  }

  /* Cluster and field masses if we were conserving mass within particles */
  const double part_clust_mass = sp->mass_init * sp->CFE;
  sp->field_mass = sp->mass_init - part_clust_mass;

  /* Determine mean cluster mass by integrating mass function,
   * along with the integral normalisation */
  double mean_mass = 0.f;
  double norm = 0.f;
  if (props->power_law_clMF) {
    /* Treat some special cases first */
    if (props->clMF_slope == 1) {
      norm = log(props->clMF_max)-log(props->clMF_min);
      mean_mass = props->clMF_max - props->clMF_min;

    } else if (props->clMF_slope == 2) {
      norm = 1./props->clMF_min - 1./props->clMF_max;
      mean_mass = log(props->clMF_max)-log(props->clMF_min);

    } else {
      norm = (pow(props->clMF_max,1.-props->clMF_slope) - 
          pow(props->clMF_min,1.-props->clMF_slope)) / 
          (1.-  props->clMF_slope);
      mean_mass = (pow(props->clMF_min,2.-props->clMF_slope) - 
          pow(props->clMF_max,2.-props->clMF_slope)) / 
          (props->clMF_slope-2.);
    }
    mean_mass /= norm;

  } else {
    /* Schechter cluster mass function */

    /* If Mcstar < clMFmin we have undefined behaviour and just form all
     * clusters at M=clMFmin */
    if (sp->Mcstar > props->clMF_min) {
      /* Get the normalisation of the integral */
      norm = schechternorm(props->clMF_min, props->clMF_max, props->clMF_slope,
          sp->Mcstar);
    }

    if (norm > 0) {
      mean_mass = meanmass(props->clMF_min, props->clMF_max, props->clMF_slope, 
          sp->Mcstar, norm);
    }
  }

  /* Number of clusters to try and form */
  if (mean_mass > 0) {
    /* Draw from Poisson distribution, so we can be resolution independent */
    sp->initial_num_clusters = poidev( (int)round(part_clust_mass/mean_mass),
        random_generator);
  }

  /* Lowest mass cluster we've formed so far */
  double min_formed_mass = FLT_MAX;
  int min_mass_idx = -1;

  /* Counter for our cluster arrays */
  int iArr = 0;

  /* draw cluster masses from mass function */
  for (int i = 0; i < sp->initial_num_clusters; i++) {
    /* Form a cluster */
    double mTry;
    if (props->power_law_clMF) {
      mTry = powerlawgen(props->clMF_min, props->clMF_max, props->clMF_slope, 
          random_generator);
    } else {
      /* Schechter cluster mass function */
      mTry = schechtergen(props->clMF_min, props->clMF_max, props->clMF_slope,
          sp->Mcstar, norm, random_generator);
    }

    sp->initial_cluster_mass_total += (float)mTry;

    /* TODO should this go within a debug statement? */
    if (!isfinite(mTry)) {
      error("mTry is not finite!! mTry=%g, Mcstar=%g, norm=%g",
          mTry, sp->Mcstar, norm);
    }

    /* Lower than our limit to evolve clusters */
    if (mTry < props->clMF_min_evolve) {
      /* Add destroyed clusters to field */
      sp->field_mass += (float)mTry;
      continue;
    }

    sp->initial_cluster_mass_evo += (float)mTry;

    /* Check if we've run out of space... */
    if (iArr < MOSAICS_MAX_CLUSTERS) {
      /* Still space left */
      //sp->clusters.id[iArr] = iArr;
      sp->clusters.mass[iArr] = (float)mTry;
      sp->clusters.initial_mass[iArr] = (float)mTry;

      /* Lowest mass cluster we've formed so far */
      if (mTry < min_formed_mass) {
        min_formed_mass = mTry;
        min_mass_idx = iArr;
      }

    } else {
      /* Ok we've run out of space. Replace the least massive
       * cluster in the array so we at least track the most massive ones */
      if (mTry > min_formed_mass) {
        /* Add the previous one to field */
        sp->field_mass += sp->clusters.mass[min_mass_idx];

        sp->clusters.mass[min_mass_idx] = (float)mTry;
        sp->clusters.initial_mass[min_mass_idx] = (float)mTry;

        /* reset the counters */
        min_formed_mass = FLT_MAX;
        min_mass_idx = -1;
        for (int j = 0; j < MOSAICS_MAX_CLUSTERS; j++) {
          if (sp->clusters.mass[j] < min_formed_mass) {
            min_formed_mass = sp->clusters.mass[j];
            min_mass_idx = j;
          }
        }

      } else {
        /* We're out of space... just add to field */
        sp->field_mass += mTry;
      }
    }
    iArr++;
  }

  /* Current number of clusters */
  sp->num_clusters = min(iArr, MOSAICS_MAX_CLUSTERS);

  /* Number of clusters above the evolution mass limit */
  sp->initial_num_clusters_evo = iArr;

  for (int i = 0; i < sp->num_clusters; i++) {
    sp->cluster_mass_total += sp->clusters.mass[i];
  }

  /* Warn when we hit the cluster array end */
  if (sp->initial_num_clusters_evo > MOSAICS_MAX_CLUSTERS) {
    message("sp->id=%lld: Trying to form more clusters (%d) than allowed (%d)\n"
        "mass_init=%g; clust_mass=%g; mean_mass=%g",
        sp->id, sp->initial_num_clusters_evo, MOSAICS_MAX_CLUSTERS, 
        sp->mass_init, part_clust_mass, mean_mass);

    /* TODO should dump some info to a file */
  }

#endif /* !defined(STAR_FORMATION_NONE) */
}

#endif /* SWIFT_MOSAICS_CLFORM_H */
