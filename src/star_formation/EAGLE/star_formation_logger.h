/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2019 Folkert Nobels (nobels@strw.leidenuniv.nl)
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
 *******************************************************************************/
#ifndef SWIFT_SCHAYE_STARFORMATION_LOGGER_H
#define SWIFT_SCHAYE_STARFORMATION_LOGGER_H

/* Some standard headers */
#include <stdlib.h>

/* Local includes */
#include "cell.h"
#include "cosmology.h"
#include "hydro.h"
#include "part.h"
#include "star_formation_logger_struct.h"

/**
 * @brief Update the star foramtion history in the current cell after creating
 * the new star particle spart sp
 *
 * @param sp new created star particle
 * @param sf the star_formation_history struct of the current cell
 */
INLINE static void star_formation_update_SFH(
    struct spart *sp, struct star_formation_history *sf) {
  /* Add mass of created sparticle to the total stellar mass in this cell*/
  sf->new_stellar_mass = sf->new_stellar_mass + sp->mass;
}

/**
 * @brief Initialize the star formation history struct
 *
 * @param sf the star_formation_history struct we want to initialize
 */
INLINE static void star_formation_init_stellar_mass(struct star_formation_history *sf) {
  /* Initialize the stellar mass to zero*/
  sf->new_stellar_mass = 0.f;
}

/**
 * @brief Initialize the star formation history struct
 *
 * @param sf the star_formation_history struct we want to initialize
 */
INLINE static void star_formation_init_SFH_active(struct star_formation_history *sf) {
  /* Initialize the active SFR */
  sf->SFR_active = 0.f;

  /* Initialize the SFR*dt active */
  sf->SFRdt_active = 0.f;

  /* Initialize the inactive SFR */
  sf->SFR_inactive = 0.f;
}

/**
 * @brief Initialize the star formation history struct
 *
 * @param sf the star_formation_history struct we want to initialize
 */
INLINE static void star_formation_init_SFH_inactive(struct star_formation_history *sf) {
  
  //message("SFR active = %e",sf->SFR_active);
  /* The active SFR becomes the inactive SFR */
  sf->SFR_inactive += sf->SFR_active;
  //message("SFR inactive = %e, SFR active = %e",sf->SFR_inactive, sf->SFR_active);

  /* Initialize the active SFR */
  sf->SFR_active = 0.f;

  /* Initialize the SFR*dt active */
  sf->SFRdt_active = 0.f;
}

/**
 * @brief function to add the progeny SFH to the parent SFH.
 *
 * @param sf parent SFH struct
 * @param sfprogeny progeny SFH struct
 */
INLINE static void star_formation_add_progeny_SFH(
    struct star_formation_history *sf,
    const struct star_formation_history *sfprogeny) {
  /* Add the new stellar mass from the progeny */
  sf->new_stellar_mass = sf->new_stellar_mass + sfprogeny->new_stellar_mass;

  /* Add the SFR of the progeny */
  sf->SFR_active += sfprogeny->SFR_active;

  /* Add the SFR * dt of the progeny */
  sf->SFRdt_active += sfprogeny->SFRdt_active;

  /* Add the Inactive SFR of the progeny */
  sf->SFR_inactive += sfprogeny->SFR_inactive;
}

/**
 * @brief Get the total star formation in this cell and add it to the star
 * formation history struct
 *
 * @param c the cell of which we want to know the star formation
 * @param sf the star formation structure to which we want to add the star
 * formation
 * @param cosmo the cosmology struct
 * @param with_cosmology if we run with cosmology
 */
INLINE static void star_formation_get_total_cell(
    struct cell *c, struct star_formation_history *sf) {
  /* Get the star formation history from the cell */
  struct star_formation_history *sfcell = &c->stars.sfh;
  sf->new_stellar_mass += sfcell->new_stellar_mass;
  
  sf->SFR_active += sfcell->SFR_active;

  sf->SFRdt_active += sfcell->SFRdt_active;

  sf->SFR_inactive += sfcell->SFR_inactive;
}

/**
 * @brief Clear the total star formation in this cell
 *
 * @param c the cell of which we want to know the star formation
 */
INLINE static void star_formation_clear_total_cell(struct cell *c) {
  /* Get the star formation history from the cell */
  //struct star_formation_history *sfcell = &c->stars.sfh;
  //sfcell->new_stellar_mass = 0.f;

  //sfcell->SFR_active = 0.f;

  //sfcell->SFRdt_active = 0.f;

  //sfcell->SFR_inactive = 0.f;
}

/**
 * @brief add the star formation to the parent cell
 *
 * @param c the cell for which we want to add the star formation
 * @param sf the combined star formation history of the progeny
 */
INLINE static void star_formation_add_to_parent_cell(
    struct cell *c, struct star_formation_history *sf) {
  /* Get the star formation history from the cell */
  struct star_formation_history *sfcell = &c->stars.sfh;
  sfcell->new_stellar_mass += sf->new_stellar_mass;

  sfcell->SFR_active += sf->SFR_active;

  sfcell->SFRdt_active += sf->SFR_active;

  sfcell->SFR_inactive += sf->SFR_inactive;
}

/**
 * @brief Initialize the star formation history structure
 *
 * @param The pointer to the star formation history structure
 * */
INLINE static void star_formation_init_SFH_engine(
    struct star_formation_history *sfh) {
  sfh->new_stellar_mass = 0.f;

  sfh->SFR_active = 0.f;

  sfh->SFRdt_active = 0.f;

  sfh->SFR_inactive = 0.f;
}

/**
 * @brief Write the final SFH to a file
 *
 * @param time the simulation time
 * @param a the scale factor
 * @param z the redshift
 * @param sf the star_formation_history struct
 */
INLINE static void star_formation_write_to_file(
    const double time, const double a, const double z,
    struct star_formation_history sf, const char* baseName) {
  static const int buffersize = 300;
  /* File name */
  char fileName[buffersize];
  snprintf(fileName, buffersize, "%s_SFH_logger.txt", baseName);

  /* Open the file */
  FILE *fp;
  fp = fopen(fileName, "a");

  /* Calculate the total SFR */
  const float totalSFR = sf.SFR_active + sf.SFR_inactive;
  fprintf(fp, "%16e %12.7f %12.7f %14e %14e %14e %14e\n", time, a, z, 
          sf.new_stellar_mass, sf.SFR_active, sf.SFRdt_active, totalSFR);
  fclose(fp);
}

/**
 * @brief Initialize the SFH logger file
 *
 * @param none
 */
INLINE static void star_formation_init_file_writer(const char* baseName) {
  /* File name */
  static const int buffersize = 300;
  char fileName[buffersize];
  snprintf(fileName, buffersize, "%s_SFH_logger.txt", baseName);

  FILE *fp;
  fp = fopen(fileName, "w");
  fprintf(fp, "# Star Formation History Logger file\n");
  fprintf(fp, "# When applicable units are in internal units\n");
  fprintf(
      fp,
      "#     Time            a            z         total M_stars  SFR (active)  SFR*dt (active)  SFR (total)\n");
  fclose(fp);
}

/**
 * @brief Add the SFR tracer to the total active SFR of this cell
 *
 * @param p the #part
 * @param xp the #xpart 
 * @param sf the SFH logger struct 
 */
INLINE static void star_formation_log_for_active_particles(
    const struct part* p, const struct xpart* xp, struct star_formation_history *sf, const double dt_star){

  /* Add the SFR to the logger file */
  sf->SFR_active += xp->sf_data.SFR;

  /* Update the active SFR*dt */
  sf->SFRdt_active += xp->sf_data.SFR * dt_star;
}

INLINE static void star_formation_log_for_inactive_particles(
    const struct part* p, const struct xpart* xp, struct star_formation_history *sf){

  /* Add the SFR to the logger file */
  sf->SFR_inactive += xp->sf_data.SFR;
}

INLINE static void star_formation_SFR_rebuilt(const struct part* p, const struct xpart* xp, struct star_formation_history *sf){
  /* Add to SFR to the sf struct */
  sf->SFR_inactive += xp->sf_data.SFR;
}

INLINE static void star_formation_recurse_SFR_rebuilt(struct cell *c, const struct cell *cp){
  struct star_formation_history *sf = &c->stars.sfh;
  const struct star_formation_history *sfp = &cp->stars.sfh;
  sf->SFR_inactive += sfp->SFR_inactive;

}

#endif /* SWIFT_SCHAYE_STARFORMATION_LOGGER_H */
