/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Folkert Nobels (nobels@strw.leidenuniv.nl)
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
#ifndef SWIFT_NONE_STAR_FORMATION_STRUCT_H
#define SWIFT_NONE_STAR_FORMATION_STRUCT_H

/* Does the model need to compute some statistics on parts
   before the start of the simulation? */
#define star_formation_compute_stats 0

/**
 * @brief Star-formation-related properties stored in the extended particle
 * data.
 */
struct star_formation_xpart_data {};

#endif /* SWIFT_NONE_STAR_FORMATION_STRUCT_H */
