/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2019 Loic Hausammann (loic.hausammann@epfl.ch)
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
#ifndef SWIFT_TASK_ORDER_EAGLE_H
#define SWIFT_TASK_ORDER_EAGLE_H

#define task_order_star_formation_before_feedback 1

/**
 * @brief Place the star formation cell at the right place in the dependency
 * graph.
 *
 * In EAGLE, star formation takes place before the feedback tasks (that are
 * launched after the stars_in task).
 *
 * @param s The #scheduler.
 * @param c The #cell on which to act.
 * @param star_resort_cell The #cell where the stars re-sorting task is in this
 * hierarchy.
 */
INLINE static void task_order_addunlock_star_formation_feedback(
    struct scheduler *s, struct cell *c, struct cell *star_resort_cell) {

  scheduler_addunlock(s, star_resort_cell->hydro.stars_resort,
                      c->stars.stars_in);
}

/**
 * @brief Place the cooling cell at the right place in the dependency
 * graph.
 *
 * The default model follows EAGLE.
 *
 * In EAGLE, the cooling takes place between the hydro and the kick2.
 *
 * @param s The #scheduler.
 * @param c The #cell on which to act.
 */
INLINE static void task_order_addunlock_cooling(struct scheduler *s,
                                                struct cell *c) {

  scheduler_addunlock(s, c->hydro.end_force, c->hydro.cooling);
  scheduler_addunlock(s, c->hydro.cooling, c->super->kick2);
}

/**
 * @brief Does a cell contain any particle that needs to do the cooling ?
 *
 * @param c The #cell.
 * @param e The #engine containing information about the current time.
 * @return 1 if the #cell contains at least an active particle, 0 otherwise.
 */
__attribute__((always_inline)) INLINE static int task_order_cell_is_active_cooling(
    const struct cell *c, const struct engine *e) {

  return cell_is_active_hydro(c, e);
}


/**
 * @brief Does this particle need to do the cooling now ?
 *
 * @param p The #part.
 * @param e The #engine containing information about the current time.
 * @return 1 if the #part is active, 0 otherwise.
 */
__attribute__((always_inline)) INLINE static int task_order_part_is_active_cooling(
    const struct part *p, const struct engine *e) {
  return part_is_active(p, e);
}

#endif /* SWIFT_TASK_ORDER_EAGLE_H */
