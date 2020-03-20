/* ***** plummer.c *****
 *
 * Create spherical Plummer model
 * Velocities are calculated using the Jeans equation for spherically symmetric
 * systems
 * Output units in M,l,t =  M_sun, pc, Myr
 *
 * Potential for the Plummer model given
 *
 *              -GM
 * Phi(r) = --------------
 *          sqrt(r^2 + b^2)
 *
 * where M is the system's total mass, and b is the Plummer scale length
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define G 4.49842e-3  // pc^3 Myr^-2 M_sun^-1
#define MASS 1e10     // M_sun
#define RADIUS 2.0    // kpc

#define TSTART 0.0

// Uncomment to check if in equilibrium
//#define CHECK_EQUILIBRIUM

int main(int argc, char *argv[]) {
  /* Generation of initial coordinates & velocities */

  double rcut = 10.0;

  int i, j;
  double ri, phi, R, ve, vr, vR, x1, x2, t1, t2, cm[6];

  FILE *file;
  int n;
  double mass, *pos, *vel;

  if (argc < 3) {
    fprintf(stderr, "Usage: %s filename nparticles\n", argv[0]);
    exit(1);
  }
  n = atoi(argv[2]);

  if (!(file = fopen(argv[1], "w"))) {
    fprintf(stderr, "Cannot open file\n");
    exit(1);
  }

  printf("file %s\n", argv[1]);
  printf("Npart %d\n", n);

  mass = MASS / (double)n;
  pos = (double *)malloc(3 * n * sizeof(double));
  vel = (double *)malloc(3 * n * sizeof(double));

  for (j = 0; j < 6; j++) cm[j] = 0.0;

  srand48(time((time_t *)NULL));

  /* Generate initial conditions from Plummer model (A & A 37, 183). */
  for (i = 0; i < n; i++) {

    /* get positions */

    do {
      x1 = drand48();
      ri = 1.0 / sqrt(pow(x1, -2.0 / 3.0) - 1.0);  // * mrpl
    } while (ri > rcut);

    phi = 2.0 * M_PI * drand48();
    pos[3 * i + 2] = (1.0 - 2.0 * drand48()) * ri;
    R = sqrt(ri * ri - pos[3 * i + 2] * pos[3 * i + 2]);
    pos[3 * i] = R * cos(phi);
    pos[3 * i + 1] = R * sin(phi);

    /* get velocities */

    /* escape velocity */
    /* ve = sqrt(2*phi) */
    ve = sqrt((2.0) / sqrt(1.0 + ri * ri));

    do {
      x2 = drand48();
      t1 = 0.1 * drand48();
      t2 = x2 * x2 * pow(1.0 - x2 * x2, 3.5);
    } while (t1 > t2);

    vr = ve * x2;

    phi = 2.0 * M_PI * drand48();
    vel[3 * i + 2] = (1.0 - 2.0 * drand48()) * vr;
    vR = sqrt(vr * vr - vel[3 * i + 2] * vel[3 * i + 2]);
    vel[3 * i] = vR * cos(phi);
    vel[3 * i + 1] = vR * sin(phi);

    /* centre of mass */
    cm[0] += pos[3 * i];
    cm[1] += pos[3 * i + 1];
    cm[2] += pos[3 * i + 2];
    cm[3] += vel[3 * i];
    cm[4] += vel[3 * i + 1];
    cm[5] += vel[3 * i + 2];
  }

  for (j = 0; j < 6; j++) cm[j] /= (double)n;

  /* scale pos and vel to analytical expectation values */
  double sx = RADIUS;  // 1/scl //3.0*M_PI/16.0;
  double sv =
      (1.02249e-6 / 1000.0) / sx *
      sqrt(pow(RADIUS * 1000.0, 3.) /
           (2.0 * 2.2489e-15 * MASS));  // (kmskpy)*scl/sct //sqrt(1.0/sx);
  for (i = 0; i < n; i++) {
    for (j = 0; j < 3; j++) {
      pos[3 * i + j] -= cm[j];
      vel[3 * i + j] -= cm[j + 3];
      pos[3 * i + j] *= sx * 1000.;    // convert to pc
      vel[3 * i + j] /= sv / 1.02273;  // convert from km/s to pc/Myr
    }
  }

  // Write particle data to file
  fprintf(file, "%d %14.7e\n", n, TSTART);
  for (i = 0; i < n; i++) {
    fprintf(file, "%d %14.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e\n", i,
            mass, pos[3 * i], pos[3 * i + 1], pos[3 * i + 2], vel[3 * i],
            vel[3 * i + 1], vel[3 * i + 2]);
  }
  fclose(file);

  double ff, tcr, rpl = RADIUS * 1000.;
  ff = 32.0 * M_PI * M_PI * sqrt(2.0) * (22.0 * 22.0 * 22.0) /
       (1.0 + 22.0 * 22.0) * rpl * rpl * rpl;
  tcr = sqrt(ff / (16.0 * G * MASS));
  printf("Crossing time at Plummer scale radius ~= %.4g Myrs\n", tcr);

#ifdef CHECK_EQUILIBRIUM
  // Initial energy
  // T = 1/2 sum_i m_i v_i^2
  // U = -1/2 sum_i,j G m_i m_j / |x_i - x_j|
  double T = 0., U = 0., x, y, z, rx, ry, rz;
  printf("Checking stability using virial theorem: 2T+U=0\n");
  for (i = 0; i < n; i++)
    T += vel[3 * i] * vel[3 * i] + vel[3 * i + 1] * vel[3 * i + 1] +
         vel[3 * i + 2] * vel[3 * i + 2];
  T *= 0.5 * mass;
  for (i = 0; i < n; i++) {
    x = pos[3 * i];
    y = pos[3 * i + 1];
    z = pos[3 * i + 2];
    for (j = i + 1; j < n; j++) {  // only count each pair once
      rx = x - pos[3 * j];
      ry = y - pos[3 * j + 1];
      rz = z - pos[3 * j + 2];
      U -= 1. / sqrt(rx * rx + ry * ry + rz * rz);
    }
  }
  U *= G * mass * mass;
  printf("T = %.4g, U = %.4g, E = %.4g, -2T/U = %.4g\n\n", T, U, T + U,
         -2. * T / U);
#endif

  return 0;
}
