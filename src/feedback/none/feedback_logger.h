/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2019 Camila Correa (correa@strw.leidenuniv.nl)
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
#ifndef SWIFT_FEEDBACK_NONE_LOGGER_H
#define SWIFT_FEEDBACK_NONE_LOGGER_H


/**
 * @brief Initialize the neutron stars merger (NSM) logger file
 *
 * @param fp the file pointer
 * @param us The current internal system of units.
 * @param phys_const Physical constants in internal units
 */
INLINE static void NSM_logger_init_log_file(
    FILE *fp, const struct unit_system *restrict us,
    const struct phys_const *phys_const) {
    
}

INLINE static void NSM_logger_write_to_log_file(
    FILE *fp, const double time,
    const struct cosmology *restrict cosmo,
    const int step, struct spart *restrict si) {
}



#endif /* SWIFT_FEEDBACK_NONE_LOGGER_H */

