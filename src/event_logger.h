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
#ifndef SWIFT_EVENT_LOGGER_H
#define SWIFT_EVENT_LOGGER_H

/* Config parameters. */
#include "../config.h"

/* Select the correct feedback model */
#if defined(FEEDBACK_NONE)
#include "./event_logger/none/event_logger.h"
#elif defined(FEEDBACK_EAGLE)
#include "./event_logger/EAGLE/event_logger.h"
#elif defined(FEEDBACK_COLIBRE)
#include "./event_logger/COLIBRE/event_logger.h"
#else
#error "Invalid choice of feedback model"
#endif

#endif /* SWIFT_EVENT_LOGGER_H */
