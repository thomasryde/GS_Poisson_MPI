/* gauss_seidel.h - Poisson problem
 *
 */
#ifndef _GAUSS_SEIDEL_H
#define _GAUSS_SEIDEL_H

// define your function prototype here
void gauss_seidel(double ***, double ***, int, int, double *);

void gauss_seidel_OMP_base(double ***, double ***, int, int, double *);

void gauss_seidel_OMP1(double ***, double ***, int, int, double *);

void gauss_seidel_OMP2(double ***, double ***, int, int, double *);

#endif
