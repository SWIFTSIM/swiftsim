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

#include <math.h>

/* Local includes */
#include "star_formation.h"
#include "physical_constants.h"
#include "dsyevj3.h"
#include "mosaics_cfelocal.h"

/**
 * @brief Do the cluster formation for the mosaics subgrid star cluster model
 *
 * @param sp The particle to act upon
 * @param stars_properties Properties of the stars model.
 * @param sf_props the star formation law properties to use
 * @param phys_const the physical constants in internal units.
 */
__attribute__((always_inline)) INLINE static void mosaics_clform(
    struct spart* restrict sp, const struct stars_props* props,
    const struct star_formation* sf_props, 
    const struct phys_const* phys_const) {

  /* TODO unit conversions into physical units */

  const double const_G = phys_const->const_newton_G;

  /* -------- Get CFE -------- */
  /* In units of kg, m, s */

  double rholoc = sp->birth_density * props->density_to_kgm3;

  double sigmaloc;
  if (props->subgrid_gas_vel_disp) {
    /* Sub-particle turbulent velocity dispersion */
    sigmaloc = sqrt(sp->birth_pressure / sp->birth_density);
  } else {
    /* The resolved velocity dispersion */
    /* 1/sqrt(3) converts to 1D assuming isotropy */
    sigmaloc = sp->gas_vel_disp / sqrt(3.f);
  }
  sigmaloc *= props->velocity_to_ms;

  /* Sound speed of cold ISM */
  double csloc;
  if (props->Fixedcs > 0) {
    csloc = props->Fixedcs;
  } else {
    csloc = sp->sound_speed_subgrid * props->velocity_to_ms;
  }

  /* Calculate CFE based on local conditions (Kruijssen 2012). units kg, m, s*/
  sp->CFE = f_cfelocal(rholoc, sigmaloc, csloc, props);

  /* Get feedback timescale while we're here and convert to internal units */
  double tfb = feedback_timescale(rholoc, sigmaloc, csloc, props);
  tfb /= props->time_to_cgs;

  /* -------- Get Mcstar -------- */

  /* Gas surface density (Krumholz & McKee 2005) */
  double phi_P = 1.f;
  if ( (sp->fgas < 1.f) && (sp->star_vel_disp > 0) ) {
    phi_P = 1.f + sp->gas_vel_disp / sp->star_vel_disp * (1.f / sp->fgas - 1.f);
  }

  const double pressure = 
      sp->birth_density * sp->gas_vel_disp * sp->gas_vel_disp / 3.f;

  const double SigmaG = sqrt(2. * pressure / (M_PI * const_G * phi_P));

  /* Toomre mass via tidal tensors (Pfeffer+18) */

  /* Temporary array for eigvec/val calculation */
  double tidevec[3][3] = {{0}}, tideval[3] = {0};
  double tide[3][3] = {{0}};
  tide[0][0] = sp->tidal_tensor[2][0];
  tide[0][1] = sp->tidal_tensor[2][1];
  tide[1][0] = sp->tidal_tensor[2][1];
  tide[0][2] = sp->tidal_tensor[2][2];
  tide[2][0] = sp->tidal_tensor[2][2];
  tide[1][1] = sp->tidal_tensor[2][3];
  tide[1][2] = sp->tidal_tensor[2][4];
  tide[2][1] = sp->tidal_tensor[2][4];
  tide[2][2] = sp->tidal_tensor[2][5];

  /* Get the eigenvectors */
  dsyevj3(tide, tidevec, tideval);

  /* sort eigenvalues, tideval[2] is largest eigenvalue */
  sort3(tideval);

  /* Circular and epicyclic frequencies */
  double Omega2 = fabs(-tideval[0] - tideval[1] - tideval[2]) / 3.;
  double kappa2 = fabs(3 * Omega2 - tideval[2]);

  sp->Omega = sqrt(Omega2);
  sp->kappa = sqrt(kappa2);

  /* 4 * pi^5 * G^2 */
  const double MTconst =
      4.f * M_PI * M_PI * M_PI * M_PI * M_PI * const_G * const_G;
  sp->Toomre_mass = MTconst * SigmaG * SigmaG * SigmaG / (kappa2 * kappa2);

  /* Toomre collapse fraction (Reina-Campos & Kruijssen 2017) */
  /* Total collapse time of Toomre volume */
  const double tcollapse = sqrt(2. * M_PI / kappa2);
  sp->frac_collapse = fmin(1., tfb / tcollapse);
  sp->frac_collapse *= sp->frac_collapse; /* f^2 */
  sp->frac_collapse *= sp->frac_collapse; /* f^4 */

  /* The 'GMC' mass */
  const double M_collapse = sp->Toomre_mass * sp->frac_collapse;

  /* Exponential truncation of cluster mass function */
  if (props->SFE > 0) {

    sp->Mcstar = props->SFE * sp->CFE * M_collapse;

  } else {

#if defined(STAR_FORMATION_COLIBRE)
    const double sfe_ff = sf_props->sfe;
#elif defined(STAR_FORMATION_GEAR)
    const double sfe_ff = sf_props->star_formation_efficiency;
#else
    /* No SFE defined for star formation model */
    /* Elmegreen 2002 */
    const double sfe_ff = 0.012;
#endif

    //TODO I'm not sure which density we want here? Same as for tfb?
    /* Free-fall time */
    const double tff = sqrt(3.0 * M_PI / (32.0 * const_G * sp->birth_density));
    //const double tff = sqrt(3.0 * M_PI / (32.0 * const_G * sp->birth_subgrid_dens));

    /* Integrated SFE. End of collapse defined by shortest of Toomre collapse
     * timescale and feedback timescale */
    double SFE_int = sfe_ff * fmin(tcollapse, tfb) / tff;

    if ((SFE_int < 0.0) || (SFE_int > 1.0)) {
      error("Integrated SFE unphysical! SFE_int=%g, sfe_ff=%g, tff=%g,"
            "tfb=%g, tcoll=%g.",SFE_int, sfe_ff, tff, tfb, tcollapse);
    }

    sp->Mcstar = SFE_int * sp->CFE * M_collapse;
  }


  /* TODO Now set up the cluster population */



  /* TODO if no clusters above mass limit were formed, set gcflag=0 */
}

#endif /* SWIFT_MOSAICS_CLFORM_H */
