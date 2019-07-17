###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2019 Loic Hausammann (loic.hausammann@epfl.ch)
#                    Josh Borrow (joshua.borrow@durham.ac.uk)
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
##############################################################################

from swiftsimio import load
from swiftsimio.visualisation import scatter
from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt

from numpy import max, min
import numpy as np

try:
    # Try and load this, otherwise we're stuck with serial
    from p_tqdm import p_map

    map = p_map
except:
    print("Try installing the p_tqdm package to make movie frames in parallel")
    pass

reverse_axis = False


def project(data, m_res, property, ylim):
    x, y, _ = data.gas.coordinates.value.T

    if reverse_axis:
        x, y = y, x

    mask = np.logical_and(y >= ylim[0], y <= ylim[1])

    x = x[mask]
    y = y[mask] - np.float64(ylim[0])

    h = data.gas.smoothing_length[mask]

    if property == "density":
        property = "masses"

    if property is not None:
        quant = getattr(data.gas, property).value[mask]
    else:
        quant = np.ones_like(x)

    image = scatter(x=x, y=y, m=quant, h=h, res=m_res)

    if reverse_axis:
        return image
    else:
        return image.T


def load_and_make_image(filename, res, property):
    image = np.zeros(res, dtype=np.float32)
    m_res = min(res)
    border = int(0.2 * m_res)

    # first part of the image
    ylim = np.array([0., 1.])

    data = load(filename)
    image[:m_res, :m_res] = project(data, m_res, property, ylim)
    if property != "density":
        image[:m_res, :m_res] /= project(data, m_res, None, ylim)

    # second part of the image
    ylim = np.array([0.5, 1.5])

    left = -m_res + border
    if reverse_axis:
        image[:, left:] = project(data, m_res, property, ylim)[:, left:]
        if property != "density":
            image[:, left:] /= project(data, m_res, None, ylim)[:, left:]
    else:
        image[left:, :] = project(data, m_res, property, ylim)[left:, :]
        if property != "density":
            image[left:, :] /= project(data, m_res, None, ylim)[left:, :]

    return image


def create_movie(filename, start, stop, resolution, property, output_filename):
    """
    Creates the movie with:

    snapshots named {filename}_{start}.hdf5 to {filename}_{stop}.hdf5
    at {resolution} x {resolution} pixel size and smoothing the given
    {property} and outputting to {output_filename}.
    """

    def baked_in_load(n):
        f = filename + "_{:04d}.hdf5".format(n)
        return load_and_make_image(f, resolution, property)

    # Make frames in parallel (reading also parallel!)
    frames = map(baked_in_load, list(range(start, stop)))

    vmax = max(list(map(max, frames)))
    vmin = min(list(map(min, frames)))

    fig, ax = plt.subplots(figsize=(1, 1))
    ax.axis("off")
    fig.subplots_adjust(0, 0, 1, 1)

    image = ax.imshow(frames[0], origin="lower", vmax=vmax, vmin=vmin)

    def frame(n):
        image.set_array(frames[n])

        return (image,)

    animation = FuncAnimation(fig, frame, range(0, stop-start), interval=40)
    animation.save(output_filename, dpi=resolution[0])


if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser(description="Creates a movie of the whole box.")

    parser.add_argument(
        "-s",
        "--snapshot",
        help="Snapshot name. Default: rayleigh_taylor",
        type=str,
        default="rayleigh_taylor",
    )

    parser.add_argument(
        "-i", "--initial", help="Initial snapshot. Default: 0",
        type=int, default=0
    )

    parser.add_argument(
        "-f", "--final", help="Final snapshot. Default: 40",
        type=int, default=40
    )

    parser.add_argument(
        "-o",
        "--output",
        help="Output filename. Default: rayleigh_taylor.mp4",
        type=str,
        default="rayleigh_taylor.mp4",
    )

    parser.add_argument(
        "-p",
        "--property",
        help="(swiftsimio) Property to plot. Default: density",
        type=str,
        default="density",
    )

    parser.add_argument(
        "-r",
        "--resolution",
        help="Resolution to make the movie at. Default: 512",
        type=int,
        default=512,
    )

    vars = parser.parse_args()

    yres = int(1.5 * vars.resolution)
    if reverse_axis:
        vars.resolution = [vars.resolution, yres]
    else:
        vars.resolution = [yres, vars.resolution]

    create_movie(
        filename=vars.snapshot,
        start=vars.initial,
        stop=vars.final,
        resolution=vars.resolution,
        property=vars.property,
        output_filename=vars.output,
    )
