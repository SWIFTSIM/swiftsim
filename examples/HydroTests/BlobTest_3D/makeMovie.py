"""
Makes a movie of the output of the blob test.

Josh Borrow (joshua.borrow@durham.ac.uk) 2019

LGPLv3
"""

from swiftsimio import load
from swiftsimio.visualisation import scatter

from p_tqdm import p_map

import matplotlib.pyplot as plt
import numpy as np

from matplotlib.colors import LogNorm
from matplotlib.animation import FuncAnimation

info_frames = 15
start_frame = 0
end_frame = 101
resolution = 1024
snapshot_name = "blob"
cmap = "Spectral_r"
text_args = dict(color="black")
# plot = "pressure"
# name = "Pressure $P$"
plot = "density"
name = "Fluid Density $\\rho$"


def get_image(n):
    """
    Gets the image for snapshot n, and also returns the associated
    SWIFT metadata object.
    """
    filename = f"{snapshot_name}_{n:04d}.hdf5"

    data = load(filename)

    x, y, _ = data.gas.coordinates.value.T

    # This is an oblong box but we want the square surrounding the cloud.

    cloud_centre_x = x[data.gas.density.argmax()]
    boxsize = data.metadata.boxsize[0].value
    half_boxsize = 0.5 * boxsize

    # Wrap the box if we need to!
    if cloud_centre_x > 0.74 * data.metadata.boxsize[0].value:
        # Boost the particles in the bottom half's x value so they wrap
        mask_lower = x <= half_boxsize
        x[mask_lower] += boxsize
        x -= half_boxsize
        cloud_centre_x -= half_boxsize
    if cloud_centre_x < 0.26 * data.metadata.boxsize[0].value:
        mask_upper = x >= half_boxsize
        x[mask_upper] -= boxsize
        x += half_boxsize
        cloud_centre_x += half_boxsize

    dx = cloud_centre_x - 0.5

    mask = np.logical_and(x > (cloud_centre_x - 0.5), x <= (cloud_centre_x + 0.5))
    x = x[mask] - dx  # To ensure this lives in [0, 1]
    y = y[mask]

    hsml = data.gas.smoothing_length.value[mask]

    if plot == "density":
        mass = data.gas.masses.value[mask]
        image = scatter(x=y, y=x, m=mass, h=hsml, res=resolution)
    else:
        quantity = getattr(data.gas, plot).value[mask]
        # Need to divide out the particle density for non-projected density quantities
        image = scatter(x=y, y=x, m=quantity, h=hsml, res=resolution) / scatter(
            x=y, y=x, m=np.ones_like(quantity), h=hsml, res=resolution
        )

    return image, data.metadata


def get_data_dump(metadata):
    """
    Gets a big data dump from the SWIFT metadata
    """

    try:
        viscosity = metadata.viscosity_info
    except:
        viscosity = "No info"

    try:
        diffusion = metadata.diffusion_info
    except:
        diffusion = "No info"

    output = (
        "$\\bf{Blob}$ $\\bf{Test}$\n\n"
        "$\\bf{SWIFT}$\n"
        + metadata.code_info
        + "\n\n"
        + "$\\bf{Compiler}$\n"
        + metadata.compiler_info
        + "\n\n"
        + "$\\bf{Hydrodynamics}$\n"
        + metadata.hydro_info
        + "\n\n"
        + "$\\bf{Viscosity}$\n"
        + viscosity
        + "\n\n"
        + "$\\bf{Diffusion}$\n"
        + diffusion
    )

    return output


def time_formatter(metadata):
    return f"$t = {metadata.t:2.2f}$"


# Generate the frames and unpack our variables
images, metadata = zip(*p_map(get_image, list(range(start_frame, end_frame))))

# The edges are funny because of the non-periodicity.
central_region = images[0][
    resolution // 10 : resolution - resolution // 10,
    resolution // 10 : resolution - resolution // 10,
]
norm = LogNorm(vmin=np.min(central_region), vmax=np.max(central_region), clip="black")

fig, ax = plt.subplots(figsize=(8, 8), dpi=resolution // 8)

fig.subplots_adjust(0, 0, 1, 1)
ax.axis("off")

# Set up the initial state
image = ax.imshow(np.zeros_like(images[0]), norm=norm, cmap=cmap, origin="lower")

description_text = ax.text(
    0.5,
    0.5,
    get_data_dump(metadata[0]),
    va="center",
    ha="center",
    **text_args,
    transform=ax.transAxes,
)

time_text = ax.text(
    0.975,
    0.975,
    time_formatter(metadata[0]),
    **text_args,
    va="top",
    ha="right",
    transform=ax.transAxes,
)

info_text = ax.text(
    0.025, 0.975, name, **text_args, va="top", ha="left", transform=ax.transAxes
)


def animate(n):
    # Display just our original frames at t < 0
    if n >= 0:
        image.set_array(images[n])
        description_text.set_text("")
        time_text.set_text(time_formatter(metadata[n]))

    return (image,)


animation = FuncAnimation(
    fig, animate, range(start_frame - info_frames, end_frame), interval=40
)

animation.save("blob.mp4")
