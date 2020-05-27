/*
* ----------------------------------------------------------------------------
* Numerical diagonalization of 3x3 matrcies
* Copyright (C) 2006  Joachim Kopp
* ----------------------------------------------------------------------------
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
* This library is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
* Lesser General Public License for more details.
*
* You should have received a copy of the GNU Lesser General Public
* License along with this library; if not, write to the Free Software
* Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
* ----------------------------------------------------------------------------


* ----------------------------------------------------------------------------
      SUBROUTINE DSYEVJ3(A, Q, W)
* ----------------------------------------------------------------------------
* Calculates the eigenvalues and normalized eigenvectors of a symmetric 3x3
* matrix A using the Jacobi algorithm.
* The upper triangular part of A is destroyed during the calculation,
* the diagonal elements are read but not destroyed, and the lower
* triangular elements are not referenced at all.
* ----------------------------------------------------------------------------
* Parameters:
*   A: The symmetric input matrix
*   Q: Storage buffer for eigenvectors
*   W: Storage buffer for eigenvalues
* ----------------------------------------------------------------------------
*/

#ifndef DSYEVJ3_H
#define DSYEVJ3_H

#include <math.h>

__attribute__((always_inline)) INLINE static void dsyevj3(double A[][3],
                                                          double W[]) {

  const int N = 3;

  /* .. Local Variables .. */
  double SD, SO;
  double S, T, U;
  double G, H, Z, THETA;
  double THRESH;
  int I, X, Y, R;

  /* Initialize Q to the identitity matrix */
  /* --- This loop is omitted since only the eigenvalues are desired --- 
  for (X = 0; X < N; X++) {
    Q[X][X] = 1.0;
    for (Y = 0; Y < X; Y++) {
      Q[X][Y] = 0.0;
      Q[Y][X] = 0.0;
    }
  }
  */

  /* Initialize W to diag(A) */
  for (X = 0; X < N; X++) W[X] = A[X][X];

  /* Calculate SQR(tr(A)) */
  SD = 0.0;
  for (X = 0; X < N; X++) SD += fabs(W[X]);
  SD *= SD;

  /* Main iteration loop */
  for (I = 0; I < 50; I++) {
    /* Test for convergence */
    SO = 0.0;
    for (X = 0; X < N; X++)
      for (Y = X + 1; Y < N; Y++) SO += fabs(A[X][Y]);
    if (SO == 0.0) return;

    if (I < 4)
      THRESH = 0.2 * SO / (N * N);
    else
      THRESH = 0.0;

    /* Do sweep */
    for (X = 0; X < N; X++) {
      for (Y = X + 1; Y < N; Y++) {
        G = 100.0 * (fabs(A[X][Y]));
        if ((I > 4) && (fabs(W[X]) + G == fabs(W[X])) &&
            (fabs(W[Y]) + G == fabs(W[Y])))
          A[X][Y] = 0.0;
        else if (fabs(A[X][Y]) > THRESH) {
          /* Calculate Jacobi transformation */
          H = W[Y] - W[X];
          if (fabs(H) + G == fabs(H))
            T = A[X][Y] / H;
          else {
            THETA = 0.5 * H / A[X][Y];
            if (THETA < 0.0)
              T = -1.0 / (sqrt(1.0 + THETA * THETA) - THETA);
            else
              T = 1.0 / (sqrt(1.0 + THETA * THETA) + THETA);
          }

          U = 1.0 / sqrt(1.0 + T * T);
          S = T * U;
          Z = T * A[X][Y];

          /* Apply Jacobi transformation */
          A[X][Y] = 0.0;
          W[X] = W[X] - Z;
          W[Y] = W[Y] + Z;
          for (R = 0; R < X; R++) {
            T = A[R][X];
            A[R][X] = U * T - S * A[R][Y];
            A[R][Y] = S * T + U * A[R][Y];
          }
          for (R = X + 1; R < Y; R++) {
            T = A[X][R];
            A[X][R] = U * T - S * A[R][Y];
            A[R][Y] = S * T + U * A[R][Y];
          }
          for (R = Y + 1; R < N; R++) {
            T = A[X][R];
            A[X][R] = U * T - S * A[Y][R];
            A[Y][R] = S * T + U * A[Y][R];
          }

          /* Update eigenvectors */
          /* --- Omitted since only the eigenvalues are desired ---
          for (R = 0; R < N; R++) {
            T = Q[R][X];
            Q[R][X] = U * T - S * Q[R][Y];
            Q[R][Y] = S * T + U * Q[R][Y];
          }
          */
        }
      }
    }
  }
  //TODO return with error
  /*!      PRINT *, "DSYEVJ3: No convergence." */
}
/* End of subroutine DSYEVJ3 */

/* Some other useful functions */

/**
 * @brief Sort array with length 3
 */
__attribute__((always_inline)) INLINE static void sort3(double val[]) {
  double tmp;
  if (val[0] > val[1]) {
    /* swap */
    tmp = val[1];
    val[1] = val[0];
    val[0] = tmp;
  }
  if (val[0] > val[2]) {
    /* swap */
    tmp = val[2];
    val[2] = val[0];
    val[0] = tmp;
  }
  if (val[1] > val[2]) {
    /* swap */
    tmp = val[2];
    val[2] = val[1];
    val[1] = tmp;
  }
}

#endif /* DSYEVJ3_H */
