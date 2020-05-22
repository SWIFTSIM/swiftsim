/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020 Evgenii Chaikin (chaikin@strw.leidenuniv.nl)
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
#ifndef SWIFT_COLIBRE_ARCLENGTH_H
#define SWIFT_COLIBRE_ARCLENGTH_H

/**
 * @brief Returns the arclength on a sphere of radius r_sphere
 * between two points with angular coordinates (theta_1, phi_1) and (theta_2,
 * phi_2)
 *
 * @param theta_1 Polar angle of point 1; \theta \in [-\pi/2, \pi/2]
 * @param phi_1 Azimuthal angle of point 1; \phi \in [-\pi, \pi)
 * @param theta_2 Polar angle of point 2; \theta \in [-\pi/2, \pi/2]
 * @param phi_2 Azimuthal angle of point 2; \phi \in [-\pi, \pi)
 * @param r_sphere Radius of the sphere on which the arclength between the two
 * points is computed
 */
__attribute__((always_inline)) INLINE static float arclength(
    const double theta_1, const double phi_1, const double theta_2,
    const double phi_2, const float r_sphere) {

  const double delta_theta = theta_2 - theta_1;
  const double delta_phi = phi_2 - phi_1;

  const double sine_theta = sin(delta_theta / 2.0);
  const double sine_phi = sin(delta_phi / 2.0);

  const float arc_length =
      asin(sqrt(sine_theta * sine_theta +
                cos(theta_1) * cos(theta_2) * sine_phi * sine_phi));

  return 2.f * r_sphere * arc_length;
}

/**
 * @brief Calculates the arclength on a sphere of radius r_sphere between
 * the ray with angular coordinates (theta_ray, phi_ray) and gas particle
 * with angular coordinates (theta_part, phi_part), and compares this calculated
 * arclength with the currently existing arclength for the considered ray
 * (the latter is known from the previous call of this function for the same ray
 * or, if the current call is the first call, is the default value)
 * Returns the new arclength if it is smaller than the older one or zero
 * otherwise
 *
 * @param theta_ray Polar angle of point 1;  \theta \in [0, \pi]
 * @param phi_ray Azimuthal angle of point 1;  \phi \in [-\pi, \pi)
 * @param theta_part Polar angle of point 2; \theta \in [0, \pi]
 * @param phi_part Azimuthal angle of point 2; \phi \in [-\pi, \pi)
 * @param r_sphere Radius of the sphere on which the arclength between the two
 * points is computed
 * @param current_archlength Current minimal value of the ray's arclength
 */
__attribute__((always_inline)) INLINE static float find_min_arclength(
    const double theta_ray, const double phi_ray, const double theta_part,
    const double phi_part, const float r_sphere,
    const float current_arclength) {

  /* We shift \theta by -\pi/2 becasue the equation we use to calculate
   * arclengths requires \theta \in [-\pi/2, \pi/2] but we have \theta \in
   * [0,\pi]  */
  const float new_arclength = arclength(
      theta_ray - M_PI_2, phi_ray, theta_part - M_PI_2, phi_part, r_sphere);

  /* Compare the arclength that has just been computed with what the ray already
   *  has. If the new one is smaller or the ray does not have any arclength yet
   * (current_arclength < 0.f), then return the new one. */
  if (new_arclength < current_arclength || current_arclength < 0.f) {
    return new_arclength;
  }
  /* If the new one is larger than the older one return zero */
  else {
    return 0.f;
  }
}

#endif /* SWIFT_COLIBRE_ARCLENGTH_H */
