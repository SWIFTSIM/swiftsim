"""
###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2017
#
# Josh Borrow (joshua.borrow@durham.ac.uk)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
###############################################################################
"""

import numpy as np
import write_gadget as wg
import h5py as h5
from scipy.special import erfinv


def inverse_gaussian(m_i, central_radius, standard_deviation):
    """
    The inverse of the Gaussian PDF used for generating the radial positions of
    the particles.

    @param: m_i | float / array-like
        - the m_i that are to be used to generate the radii of the particles

    @param: central_radius | float
        - central radius of the guassian

    @param: standard_deviation | float
        - standard deviation of the gaussian.

    ---------------------------------------------------------------------------

    @return: radius | float / array-like
        - the radius associated with m_i (see README for theory).
    """
    error_function = erfinv(2 * m_i - 1)
    radius = central_radius + standard_deviation * np.sqrt(2) * error_function

    return radius


def generate_m_i(n_particles):
    """
    Generate the m_i for each particle

    @param: n_particles | int
        - number of m_i to generate, or equivalently the number of particles.

    ---------------------------------------------------------------------------

    @return: m_i | numpy.array
        - the m_i that are used to generate the radii of each particle.
    """
    m_i = (np.arange(n_particles) + 0.5) / n_particles

    return m_i


def generate_theta_i(r_i, theta_initial=0.):
    """
    Generate the theta associated with the particles based on their radii.

    @param: r_i | float / array-like
        - the radii of the particles.

    @param: theta_initial | float | optional
        - initial radius to start the generation at.

    ---------------------------------------------------------------------------

    @return: theta_i | numpy.array
        - angles associated with the particles in the plane.
    """
    radii_fraction = r_i[:-1] / r_i[1:]
    d_theta_i = np.sqrt(2 * np.pi * (1 - radii_fraction))

    theta_i = np.empty_like(r_i)
    theta_i[0] = theta_initial

    # Need to do this awful loop because each entry relies on the one before it.
    # Unless someone knows how to do this in numpy.
    for i in range(len(d_theta_i)):  # first is theta_initial
        theta_i[i+1] = theta_i[i] + d_theta_i[i]

    return theta_i


def convert_polar_to_cartesian(r_i, theta_i):
    """
    Calculate the cartesian co-ordinates (to be used to store in the GADGET
    file) from the polar co-ordinates generated by the script.abs

    @param: r_i | float / array-like
        - the radii of the particles

    @param: theta_i | float / array-like
        - the polar angle co-ordinate of the particles

    ---------------------------------------------------------------------------

    @return: x_i | float / array-like
        - the x co-ordinates of the particles

    @return: y_i | float / array-like
        - the y co-ordinates of the particles
    """
    x_i = r_i * np.cos(theta_i)
    y_i = r_i * np.sin(theta_i)

    return x_i, y_i


def get_keplerian_velocity(r_i, theta_i, mass):
    r"""
    The keplerian velocity of the particles given their position in the ring.
    Note that this is exactly in the direction of $\hat{theta}$.

    @param: r_i | float / array-like
        - radial positions of the particles

    @param: theta_i | float / array_like
        - polar angle of the particles

    @param: mass | float
        - GM (i.e. Newton's Gravitational Constant multiplied by the mass) of
          the central point particle that the keplerian ring rotates around.

    ---------------------------------------------------------------------------

    @return: v_x_i | numpy.array
        - the velocity of particles in the x direction

    @return: v_y_i | numpy.array
        - the velocity of particles in the y direction.
    """
    force_modifier = np.sqrt(mass / r_i)

    v_x_i = force_modifier * (-np.sin(theta_i))
    v_y_i = force_modifier * (np.cos(theta_i))

    return v_x_i, v_y_i


def QSP_fix(r_i, theta_i):
    """
    The start and end of the disk will have the end of the spiral there. That
    is no good and introduces some shear forces, so we need to move them to
    concentric circles. We'll also add some extra particles on the final
    'layer' of this giant onion to ensure that all of the circles are complete.

    @param: r_i | numpy.array
        - the original r_i generated from inverse_gaussian (and perhaps
          afterwards masked).

    @param: theta_i | numpy.array
        - the original theta_i generated from generate_theta_i.

    ---------------------------------------------------------------------------

    @return r_i_fixed | numpy.array
        - the fixed, concentric circle-like r_i. Note that these arrays do not
          necessarily have the same length as the r_i, theta_i that are
          passed in and you will have to re-calculate the number of particles
          in the system.

    @return theta_i_fixed | numpy.array
        - the fixed, concentric circle-like theta_i
    """

    # Thankfully, the delta_thetas are not wrapped (i.e. they keep on going
    # from 2 pi -- we need to split the arrays into 'circles' by splitting
    # the theta_i array every 2 pi.

    rotations = 1
    circles = []

    these_r_i = []
    these_theta_i = []

    for radius, theta in zip(r_i, theta_i):
        if theta > rotations * 2 * np.pi:
            circles.append([these_r_i, these_theta_i])
            these_r_i = []
            these_theta_i = []
            rotations += 1

        these_r_i.append(radius)
        these_theta_i.append(theta)

    # Now we need to do the averaging.
    # We want to have all particles in each circle a fixed radius away from the
    # centre, as well as having even spacing between each particle. The final
    # ring may be a bit dodgy still, but we will see.
    
    r_i_fixed = []
    theta_i_fixed = []

    for circle in circles:
        n_particles = len(circle[0])
        radius = sum(circle[0]) / n_particles
        theta_initial = circle[1][0]

        theta_sep = 2 * np.pi / n_particles
        
        theta = [t * theta_sep for t in range(n_particles)]
        radii = [radius] * n_particles

        r_i_fixed += radii
        theta_i_fixed += theta

    return np.array(r_i_fixed), np.array(theta_i_fixed)


def generate_particles(n_particles, central_radius, std_dev, mass, int_e):
    """
    A quick wrapper function that generates the x and y co-ordinates of
    particles in a keplerian ring.

    @param: n_particles | int
        - the number of particles in the ring

    @param: central_radius | float
        - the radius around which the particles are arranged

    @param: std_dev | float
        - the standard deviation of the gaussian which determines how the
          particles are arranged horizontally and vertically around the ring.

    @param: mass | float
        - GM (i.e. Newton's Gravitational Constant multiplied by the mass) of
          the central point particle that the keplerian ring rotates around.

    @param: int_e | float
		- the internal energy of each individual particle.

    ---------------------------------------------------------------------------

    @return: x_i | numpy.array
        - the x co-ordinates of the particles in the ring

    @return: y_i | numpy.array
        - the y co-ordinates of the particles in the ring

    @return: v_x_i | numpy.array
        - the velocity of particles in the x direction

    @return: v_y_i | numpy.array
        - the velocity of particles in the y direction.
    """
    m_i = generate_m_i(n_particles)
    r_i = inverse_gaussian(m_i, central_radius, std_dev)
    theta_i = generate_theta_i(r_i)

    # We need to remove the start and end of the spiral. We can do that here.
    r_i, theta_i = QSP_fix(r_i, theta_i)

    # We need to remove those that are too large or small and end up as noisy
    # neighbors and really disrupt the ring.

    upper_bound = r_i > central_radius + 1.5 * std_dev
    lower_bound = r_i < central_radius - 1.5 * std_dev

    mask = np.logical_or(upper_bound, lower_bound)

    # Masked arrays.compressed() simply returns the non-masked data.
    r_i = np.ma.masked_array(r_i, mask).compressed()
    theta_i = np.ma.masked_array(theta_i, mask).compressed()

    # Now we can continue again!

    x_i, y_i = convert_polar_to_cartesian(r_i, theta_i)
    v_x_i, v_y_i = get_keplerian_velocity(r_i, theta_i, mass)

    n_particles = len(r_i)
    int_energy = np.zeros(n_particles) + int_e

    return x_i, y_i, v_x_i, v_y_i, int_energy


def save_to_gadget(filename, x_i, y_i, v_x_i, v_y_i, int_energy, pm, hsml):
    """ Save the particle data to a GADGET .hdf5 file.

    @param: filename | string
        - filename of the hdf5 file to save.

    @param: x_i | array-like
        - x positions of the particles

    @param: y_i | array-like
        - y positions of the particles

    @param: v_x_i | array-like
        - x velocities of the particles

    @param: v_y_i | array-like
        - y velocities of the particles

    @param: hsml | float
        - smoothing length of the particles.

    @param: pm | float
        - mass of the particles.

    ---------------------------------------------------------------------------
    """
    n_particles = len(x_i)

    positions = np.array([x_i, y_i, np.zeros_like(x_i)]).T
    velocities = np.array([v_x_i, v_y_i, np.zeros_like(v_x_i)]).T
    ids = np.arange(n_particles)

    with h5.File(filename, "w") as handle:
        wg.write_header(
            handle,
            boxsize=200.,
            flag_entropy=0,
            np_total=np.array([n_particles, 0, 0, 0, 0, 0]),
            np_total_hw=np.array([0, 0, 0, 0, 0, 0]),
            other={"MassTable" : np.array([pm, 0, 0, 0, 0, 0])}
        )

        wg.write_runtime_pars(
            handle,
            periodic_boundary=1,
        )

        wg.write_units(
            handle,
            current=1.,
            length=3.0856776e21,
            mass=1.9885e33,
            temperature=1.,
            time=3.085678e16,
        )

        wg.write_block(
            handle,
            0,  # gas
            positions,
            velocities,
            ids,
            mass=np.zeros(n_particles) + pm,
            int_energy=int_energy,
            smoothing=np.zeros(n_particles) + hsml
        )

    return


if __name__ == "__main__":
    # Run in script mode!

    import argparse as ap

    PARSER = ap.ArgumentParser(
        description="""Initial conditions generator for the Keplerian Ring
                       example. It has sensible defaults, but if you wish to
                       play around this script can be configured using command
                       line arguments. For more information on those, see
                       below."""
    )

    PARSER.add_argument(
        '-m',
        '--gravity_mass',
        help='GM for your central point mass (default: 43009.7 simulation units)',
        required=False,
        default=43009.7
    )
    PARSER.add_argument(
        '-f',
        '--filename',
        help='Output filename (default: initial_conditions.hdf5)',
        required=False,
        default="initial_conditions.hdf5"
    )
    PARSER.add_argument(
        '-n',
        '--nparts',
        help='Total number of particles (default: 10000). This is approximate.',
        required=False,
        default=10000
    )
    PARSER.add_argument(
        '-r',
        '--centralradius',
        help='Distance to centre of the gaussian from (0, 0) \
             (i.e. the mean radius of the ring). \
             (default: 10 simulation units)',
        required=False,
        default=10
    )
    PARSER.add_argument(
        '-sd',
        '--standarddev',
        help='Standard deviation of the gaussian (i.e. the width of the\
              gaussian) (default: 2.5 simulation units).',
        required=False,
        default=2.5
    )
    PARSER.add_argument(
        '-pm',
        '--particlemass',
        help='Mass of the SPH particles (default: 1 simulation units).',
        required=False,
        default=1
    )
    PARSER.add_argument(
        '-hs',
        '--smoothing',
        help='Smoothing length of the SPH particles (default: 1.28 simulation units).',
        required=False,
        default=1.28
    )
    PARSER.add_argument(
        '-ie',
        '--internalenergy',
        help='Internal energy of the SPH particles (default: 0.015 simulation units).',
        required=False,
        default=0.015
    )

    ARGS = vars(PARSER.parse_args())

    PARTICLES = generate_particles(
        int(ARGS['nparts']),
        float(ARGS['centralradius']),
        float(ARGS['standarddev']),
        float(ARGS['gravity_mass']),
		float(ARGS['internalenergy'])
    )

    save_to_gadget(
        ARGS['filename'],
        PARTICLES[0],
        PARTICLES[1],
        PARTICLES[2],
        PARTICLES[3],
        PARTICLES[4],
        float(ARGS['particlemass']),
        float(ARGS['smoothing'])
    )
