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
#include "cosmology.h"
#include "dsyevj3.h"
#include "engine.h"
#include "mosaics_cfelocal.h"

/**
 * @brief Do the cluster formation for the mosaics subgrid star cluster model
 *
 * @param sp The particle to act upon
 * @param stars_properties Properties of the stars model.
 * @param e The #engine
 * @param cosmo the cosmological parameters and properties.
 * @param with_cosmology if we run with cosmology.
 */
__attribute__((always_inline)) INLINE static void mosaics_clform(
    struct spart* restrict sp, const struct stars_props* props,
    const struct engine* e, const struct cosmology* cosmo,
    const int with_cosmology) {

  const double const_G = e->physical_constants->const_newton_G;
  const double const_boltzmann_k = e->physical_constants->const_boltzmann_k;
  const double const_proton_mass = e->physical_constants->const_proton_mass;

  /* molecular weight of primordial gas (mu=1.22)*/
  const double mu_neutral = e->hydro_properties->mu_neutral;

  /* TODO unit conversions into physical units */

  /* -------- Get CFE -------- */

  /* Sub-particle turbulent velocity dispersion */
  /* sqrt(3) to convert to 3D */
  double turbVelDisp = sqrt(3.f * sp->birth_pressure / sp->birth_density);

  /* TODO make this optional */
  /* Combined resolved and unresolved velocity dispersion */
  // totalVelDisp = sqrt( sp->gasVelDisp * sp->gasVelDisp +
  //                      turbVelDisp * turbVelDisp );
  double totalVelDisp = turbVelDisp;

  /* In units of kg, m, s */
  /* 1/sqrt(3) converts to 1D assuming isotropy */
  double sigmaloc = totalVelDisp / sqrt(3.f) * props->velocity_to_ms;
  double rholoc = sp->birth_density * props->density_to_kgm3;

  // TODO We might want the subgrid T for this?
  // double csloc = props->Fixedcs;
  double csloc = sp->birth_subgrid_temp * const_boltzmann_k /
                 (mu_neutral * const_proton_mass);
  csloc = sqrt(csloc) * props->velocity_to_ms;

  /* Calculate CFE based on local conditions (Kruijssen 2012). units kg, m, s*/
  sp->CFE = f_cfelocal(rholoc, sigmaloc, csloc, props);

  /* Get feedback timescale while we're here */
  double tfb = feedback_timescale(rholoc, sigmaloc, csloc, props);

  /* -------- Get Mcstar -------- */

  /* Gas surface density (Krumholz & McKee 2005) */
  double phi_P = 1.f;
  if (sp->starVelDisp > 0) {
    phi_P = 1.f + sp->gasVelDisp / sp->starVelDisp * (1.f / sp->fgas - 1.f);
  }

  double SigmaG = sqrt(2. * sp->birth_pressure / (M_PI * const_G * phi_P));

  /* Toomre mass via tidal tensors (Pfeffer+18) */

  /* Temporary array for eigvec/val calculation */
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

  double tidevec[3][3] = {{0}}, tideval[3] = {0};
  dsyevj3(tide, tidevec, tideval);
  /* sort eigenvalues, tideval[2] is largest eigenvalue */
  sort3(tideval);

  double Omega2 = fabs(-tideval[0] - tideval[1] - tideval[2]) / 3.;
  double kappa2 = fabs(3 * Omega2 - tideval[2]);

  sp->Omega = sqrt(Omega2);
  sp->kappa = sqrt(kappa2);

  /* 4 * pi^5 * G^2 */
  const double MTconst =
      4.f * M_PI * M_PI * M_PI * M_PI * M_PI * const_G * const_G;
  sp->Toomre_mass = MTconst * SigmaG * SigmaG * SigmaG / (kappa2 * kappa2);

  /* Toomre collapse fraction (Reina-Campos & Kruijssen 2017) */
  double tff = sqrt(2. * M_PI / kappa2);
  sp->fracCollapse = fmin(1., tfb / tff);
  sp->fracCollapse *= sp->fracCollapse; /* f^2 */
  sp->fracCollapse *= sp->fracCollapse; /* f^4 */

  /* The 'GMC' mass */
  double M_collapse = sp->Toomre_mass * sp->fracCollapse;

  /* Exponential truncation of cluster mass function */
  /* TODO here we should use SFE = starform->star_formation_efficiency if it
   * exists */
  sp->Mcstar = props->SFE * sp->CFE * M_collapse;

  /* TODO if no clusters above mass limit were formed, set gcflag=0 */
}

#endif /* SWIFT_MOSAICS_CLFORM_H */
