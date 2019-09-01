/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2019 Folkert Nobels (nobels@strw.leidenuniv.nl)
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
#ifndef SWIFT_NONE_FEEDBACK_LOGGER_H
#define SWIFT_NONE_FEEDBACK_LOGGER_H

/**
 * @brief Initialize the SFH logger file
 *
 * @param fp the file pointer
 * @param us The current internal system of units.
 * @param phys_const Physical constants in internal units
 */
INLINE static void feedback_SNIa_logger_init_log_file(
    FILE *fp, const struct unit_system *restrict us,
    const struct phys_const *phys_const) {}

INLINE static void feedback_SNIa_logger_write_to_log_file(
    FILE *fp, const double time, struct spart *restrict si,
    struct part *restrict pj, struct xpart *restrict xpj,
    const struct cosmology *restrict cosmo, const int step) {}

#endif /* SWIFT_NONE_FEEDBACK_LOGGER_H */
