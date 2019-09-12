#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "chimes_vars.h"
#include "chimes_proto.h"

/*
 * ----------------------------------------------------------------------
 * The following routine takes a value, x, and finds its position, i, 
 * in a 1D table, along with the discplacement, dx, of x in between 
 * the discrete table values. We assume here that the table is 
 * evenly spaced. 
 * ----------------------------------------------------------------------
 */

void chimes_get_table_index(ChimesFloat *table, int ntable, ChimesFloat x, int *i, ChimesFloat *dx)
{
  ChimesFloat denominator;

  denominator = (table[ntable - 1] - table[0]) / (ntable - 1.0);

  if(x <= table[0])
    {
      *i = 0;
      *dx = 0.0;
    }
  else if(x >= table[ntable - 1])
    {
      *i = ntable - 2;
      *dx = 1.0;
    }
  else
    {
      *i = (int) floor((x - table[0]) / denominator);
      *dx = (x - table[*i]) / denominator;
    }
}

/*
 * ----------------------------------------------------------------------
 * Routine to perform a linear interpolation
 * ----------------------------------------------------------------------
 */

ChimesFloat chimes_interpol_1d(ChimesFloat *table, int i, ChimesFloat dx)
{
  return (1 - dx) * table[i] + dx * table[i + 1];
}

/*
 * ----------------------------------------------------------------------
 * Routine to perform a bi-linear interpolation
 * ----------------------------------------------------------------------
 */

ChimesFloat chimes_interpol_2d(ChimesFloat **table, int i, int j, ChimesFloat dx, ChimesFloat dy)
{
  ChimesFloat output, dx_m, dy_m;

  dx_m = 1.0 - dx;
  dy_m = 1.0 - dy; 

  output =
    dx_m * dy_m * table[i][j] +
    dx_m * dy * table[i][j + 1] +
    dx * dy_m * table[i + 1][j] +
    dx * dy * table[i + 1][j + 1];

  return output;
}

/*
 * ----------------------------------------------------------------------
 * Routine to perform a tri-linear interpolation
 * ----------------------------------------------------------------------
 */
ChimesFloat chimes_interpol_3d(ChimesFloat ***table, int i, int j, int k, ChimesFloat dx, ChimesFloat dy, ChimesFloat dz)
{
  ChimesFloat output, dx_m, dy_m, dz_m;

  dx_m = 1.0 - dx;
  dy_m = 1.0 - dy; 
  dz_m = 1.0 - dz; 

  output =
    dx_m * dy_m * dz_m * table[i][j][k] +
    dx_m * dy_m * dz * table[i][j][k + 1] +
    dx_m * dy * dz_m * table[i][j + 1][k] +
    dx_m * dy * dz * table[i][j + 1][k + 1] +
    dx * dy_m * dz_m * table[i + 1][j][k] +
    dx * dy_m * dz * table[i + 1][j][k + 1] +
    dx * dy * dz_m * table[i + 1][j + 1][k] +
    dx * dy * dz * table[i + 1][j + 1][k + 1];

  return output;
}



/*
 * ----------------------------------------------------------------------
 * Routine to perform a quadri-linear interpolation
 * ----------------------------------------------------------------------
 */
ChimesFloat chimes_interpol_4d(ChimesFloat ****table, int i, int j, int k, int l, ChimesFloat dx, ChimesFloat dy, ChimesFloat dz, ChimesFloat dw)
{
  double output, dx_m, dy_m, dz_m, dw_m;

  dx_m = 1.0 - dx;
  dy_m = 1.0 - dy; 
  dz_m = 1.0 - dz; 
  dw_m = 1.0 - dw; 
  
  output =
    dx_m * dy_m * dz_m * dw_m * table[i][j][k][l] +
    dx_m * dy_m * dz_m * dw * table[i][j][k][l+1] +
    dx_m * dy_m * dz * dw_m * table[i][j][k+1][l] +
    dx_m * dy_m * dz * dw * table[i][j][k+1][l+1] +
    dx_m * dy * dz_m * dw_m * table[i][j+1][k][l] +
    dx_m * dy * dz_m * dw * table[i][j+1][k][l+1] +
    dx_m * dy * dz * dw_m * table[i][j+1][k+1][l] +
    dx_m * dy * dz * dw * table[i][j+1][k+1][l+1] +
    dx * dy_m * dz_m * dw_m * table[i+1][j][k][l] +
    dx * dy_m * dz_m * dw * table[i+1][j][k][l+1] +
    dx * dy_m * dz * dw_m * table[i+1][j][k+1][l] +
    dx * dy_m * dz * dw * table[i+1][j][k+1][l+1] +
    dx * dy * dz_m * dw_m * table[i+1][j+1][k][l] +
    dx * dy * dz_m * dw * table[i+1][j+1][k][l+1] +
    dx * dy * dz * dw_m * table[i+1][j+1][k+1][l] + 
    dx * dy * dz * dw * table[i+1][j+1][k+1][l+1];

  return output;
}
