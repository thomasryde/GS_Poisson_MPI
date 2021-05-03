
#include "Gauss_Seidel_Blocked.c"
void Gauss_Seidel_Blocked(double ***f,double *** u, int n, int N,int max_iter,double * tolerance);

#include "Gauss_Seidel_nonblocked.c"
void Gauss_Seidel_nonblocked(double ***f,double *** u, int n, int N,int max_iter,double * tolerance);

#include "Gauss_seidel_redblack.c"
void Gauss_seidel_redblack(double ***f,double *** u, int n, int N,int max_iter,double * tolerance);

#include "Gauss_seidel_redblack_timing.c"
void Gauss_seidel_redblack_timing(double ***f,double *** u, int n, int N,int max_iter,double * tolerance);

#include "Gauss_seidel_redblack_mp.c"
void Gauss_seidel_redblack_mp(double ***f,double *** u, int n, int N,int max_iter,double * tolerance);