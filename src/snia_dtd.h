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
#ifndef SWIFT_SNIA_DTD_H
#define SWIFT_SNIA_DTD_H

/**
 * @file src/snia_dtd.h
 * @brief Branches between the different SNIa delay time distributions recipies.
 */

/* Config parameters. */
#include "../config.h"

/* Import the right star formation law definition */
#if defined(SNIA_DTD_NONE)
#include "./feedback/DTD/none/dtd.h"
#elif defined(SNIA_DTD_EXP)
#include "./feedback/DTD/exponential/dtd.h"
#elif defined(SNIA_DTD_POWER)
#include "./feedback/DTD/powerlaw/dtd.h"
#elif defined(SNIA_DTD_POWER_BETA1)
#include "./feedback/DTD/powerlaw1/dtd.h"
#elif defined(SNIA_DTD_GAUSSIAN)
#include "./feedback/DTD/gaussian/dtd.h"
#elif defined(SNIA_DTD_CONST)
#include "./feedback/DTD/constant/dtd.h"
#else
#error "Invalid choice of star formation law"
#endif

#endif /* SWIFT_SNIA_DTD_H */
