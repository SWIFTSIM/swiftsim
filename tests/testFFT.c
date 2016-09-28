/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (C) 2015 Matthieu Schaller (matthieu.schaller@durham.ac.uk).
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

#include "../config.h"

/* Some standard headers. */
#include <stdlib.h>
#include <string.h>

#ifndef HAVE_FFTW

int main() { return 0; }

#else

#include <fftw3.h>

/* Includes. */
#include "swift.h"

const double G = 1.;

const size_t N = 16;
const size_t PMGRID = 8;

// const double asmth = 2. * M_PI * const_gravity_a_smooth / boxSize;
// const double asmth2 = asmth * asmth;
// const double fact = G / (M_PI * boxSize) * (1. / (2. * boxSize / PMGRID));

void printGrid(double* grid, const char* fname, size_t dimx, size_t dimy,
               size_t dimz) {

  FILE* file = fopen(fname, "w");
  for (size_t i = 0; i < dimx; ++i) {
    fprintf(file, "--- %zu ---\n", i);
    for (size_t j = 0; j < dimy; ++j) {
      for (size_t k = 0; k < dimz; ++k) {
        fprintf(file, "%f ", grid[i * dimx * dimy + j * dimx + k]);
      }
      fprintf(file, "\n");
    }
    fprintf(file, "\n\n");
  }
  fclose(file);
}

void printGrid_r(fftw_complex* grid, const char* fname, size_t dimx,
                 size_t dimy, size_t dimz) {

  FILE* file = fopen(fname, "w");
  for (size_t i = 0; i < dimx; ++i) {
    fprintf(file, "--- %zu ---\n", i);
    for (size_t j = 0; j < dimy; ++j) {
      for (size_t k = 0; k < dimz; ++k) {
        fprintf(file, "%e ", grid[i * dimx * dimy + j * dimx + k][0]);
      }
      fprintf(file, "\n");
    }
    fprintf(file, "\n\n");
  }
  fclose(file);
}

int main() {

  /* Initialize CPU frequency, this also starts time. */
  unsigned long long cpufreq = 0;
  clocks_set_cpufreq(cpufreq);

  /* Simulation properties */
  const size_t count = N * N * N;
  const double boxSize = 1.;

  /* Create some particles */
  struct gpart* gparts = malloc(count * sizeof(struct gpart));
  bzero(gparts, count * sizeof(struct gpart));
  for (size_t i = 0; i < N; ++i) {
    for (size_t j = 0; j < N; ++j) {
      for (size_t k = 0; k < N; ++k) {

        struct gpart* gp = &gparts[i * N * N + j * N + k];

        gp->x[0] = i * boxSize / N + boxSize / (2 * N);
        gp->x[1] = j * boxSize / N + boxSize / (2 * N);
        gp->x[2] = k * boxSize / N + boxSize / (2 * N);

        gp->mass = k * 1. / count;

        gp->id_or_neg_offset = i * N * N + j * N + k;
      }
    }
  }

  /* Properties of the mesh */
  const size_t meshmin[3] = {0, 0, 0};
  const size_t meshmax[3] = {PMGRID - 1, PMGRID - 1, PMGRID - 1};

  const size_t dimx = meshmax[0] - meshmin[0] + 2;
  const size_t dimy = meshmax[1] - meshmin[1] + 2;
  const size_t dimz = meshmax[2] - meshmin[2] + 2;

  const double fac = PMGRID / boxSize;
  // const size_t PMGRID2 = 2 * (PMGRID / 2 + 1);

  message("dimx=%zd dimy=%zd dimz=%zd", dimx, dimy, dimz);

  /* Allocate and empty the rhogrid mesh */
  const size_t rhogrid_size = dimx * dimy * dimz;
  double* rhogrid = fftw_malloc(rhogrid_size * sizeof(double));
  bzero(rhogrid, rhogrid_size * sizeof(double));

  /* Do CIC with the particles */
  for (size_t pid = 0; pid < count; ++pid) {

    const struct gpart* const gp = &gparts[pid];

    const size_t slab_x =
        (fac * gp->x[0] >= PMGRID) ? PMGRID - 1 : fac * gp->x[0];
    const size_t slab_y =
        (fac * gp->x[1] >= PMGRID) ? PMGRID - 1 : fac * gp->x[1];
    const size_t slab_z =
        (fac * gp->x[2] >= PMGRID) ? PMGRID - 1 : fac * gp->x[2];

    const double dx = fac * gp->x[0] - (double)slab_x;
    const double dy = fac * gp->x[1] - (double)slab_y;
    const double dz = fac * gp->x[2] - (double)slab_z;

    const size_t slab_xx = slab_x + 1;
    const size_t slab_yy = slab_y + 1;
    const size_t slab_zz = slab_z + 1;

    rhogrid[(slab_x * dimy + slab_y) * dimz + slab_z] +=
        gp->mass * (1.0 - dx) * (1.0 - dy) * (1.0 - dz);
    rhogrid[(slab_x * dimy + slab_yy) * dimz + slab_z] +=
        gp->mass * (1.0 - dx) * dy * (1.0 - dz);
    rhogrid[(slab_x * dimy + slab_y) * dimz + slab_zz] +=
        gp->mass * (1.0 - dx) * (1.0 - dy) * dz;
    rhogrid[(slab_x * dimy + slab_yy) * dimz + slab_zz] +=
        gp->mass * (1.0 - dx) * dy * dz;
    rhogrid[(slab_xx * dimy + slab_y) * dimz + slab_z] +=
        gp->mass * (dx) * (1.0 - dy) * (1.0 - dz);
    rhogrid[(slab_xx * dimy + slab_yy) * dimz + slab_z] +=
        gp->mass * (dx)*dy * (1.0 - dz);
    rhogrid[(slab_xx * dimy + slab_y) * dimz + slab_zz] +=
        gp->mass * (dx) * (1.0 - dy) * dz;
    rhogrid[(slab_xx * dimy + slab_yy) * dimz + slab_zz] +=
        gp->mass * (dx)*dy * dz;
  }

  for (size_t i = 0; i < dimx * dimy; ++i) rhogrid[i] = 100.f;

  /* Show the grid */
  printGrid(rhogrid, "rho.dat", dimx, dimy, dimz);

  /* Prepare the force grid */
  const size_t rhotilde_size = rhogrid_size;
  fftw_complex* rhotilde = fftw_malloc(rhotilde_size * sizeof(fftw_complex));
  bzero(rhotilde, rhotilde_size * sizeof(fftw_complex));
  for (size_t i = 0; i < dimx * dimy * dimz; ++i) rhotilde[i][0] = 1.f;

  /* Prepare the FFT */
  fftw_plan p =
      fftw_plan_dft_r2c_3d(dimx, dimy, dimz, rhogrid, rhotilde, FFTW_MEASURE);

  fftw_print_plan(p);
  printf("\n");

  /* Execute the FFT */
  fftw_execute(p);

  /* Show the FFT */
  printGrid_r(rhotilde, "rhotilde_r.dat", dimx, dimy, dimz);
  // printGrid(rhogrid, "rhotilde_r.dat", dimx, dimy, dimz);

  /* Clean-up */
  fftw_destroy_plan(p);
  fftw_free(rhotilde);
  fftw_free(rhogrid);
  free(gparts);
  return 0;
}

#endif
