/* main.c - Poisson problem in 3D
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include "alloc3d.h"
#include "print.h"
#include <math.h>
#include "InitVari.h"
#include <time.h>
#include <omp.h>

#ifdef _JACOBI
#include "jacobi.h"
#endif

#ifdef _GAUSS_SEIDEL
#include "gauss_seidel.h"
#endif

#define N_DEFAULT 100

#define TWK_OMP2 1
// #define TWK_OMP_GS_BASE 1 

int main(int argc, char *argv[]) {

    int 	N = N_DEFAULT;
    int 	iter_max = 1000;
    double 	tolerance;
    double	start_T;
    int		output_type = 0;
    char	*output_prefix = "poisson_res";
    char    *output_ext    = "";
    char	output_filename[FILENAME_MAX];
    double 	***u = NULL;
    double 	***u_old = NULL;
    double 	***f = NULL;

    /* get the paramters from the command line */
    N         = atoi(argv[1]);	// grid size
    iter_max  = atoi(argv[2]);  // max. no. of iterations
    tolerance = atof(argv[3]);  // tolerance
    start_T   = atof(argv[4]);  // start T for all inner grid points
    if (argc == 6) {
	output_type = atoi(argv[5]);  // ouput type
    }

    // allocate memory
    if ( (u = d_malloc_3d(N+2, N+2, N+2)) == NULL ) {
        perror("array u: allocation failed");
        exit(-1);
    }
    
    if ( (f = d_malloc_3d(N+2, N+2, N+2)) == NULL ) {
        perror("array f: allocation failed");
        exit(-1);
    }

    if ( (u_old = d_malloc_3d(N+2, N+2, N+2)) == NULL ) {
        perror("array u_old: allocation failed");
        exit(-1);
    }
    
    double ts, te;

    #ifdef TWK_OMP1
        printf("hello BASELINE!\n");
        InitU(u, N);
        InitU(u_old, N);
        InitF(f, N);

        ts = omp_get_wtime();        
        jacobi_OMP1(f,u,u_old,N,iter_max,&tolerance);
        te = omp_get_wtime() - ts;

    #elif TWK_OMP2
        printf("hello OMP2!\n");
        
        Init_OMP2(u, u_old, f, N);
        
        ts = omp_get_wtime();
            
        jacobi_OMP2(f,u,u_old,N,iter_max,&tolerance);      
            
        te = omp_get_wtime() - ts;
    
    
    #elif TWK_OMP_GS_BASE
        printf("hello OMP GS Baseline!\n");
        
        
        Init_OMP2(u, u_old, f, N);
        
        ts = omp_get_wtime();
            
        gauss_seidel_OMP_base(f,u,N,iter_max,&tolerance);      
            
        te = omp_get_wtime() - ts;

        int t1 = 1;
        int t2 = 1;
        int t3 = (int) ceil((double)(N+2)/3);

        double value = -6*u[t1][t2][t3] + u[t1][t2][t3 + 1] + u[t1][t2][t3-1] + u[t1][t2-1][t3] + u[t1][t2+1][t3] + u[t1-1][t2][t3] + u[t1+1][t2][t3];
        value *= (N+1)*(N+1)/4;
        value += f[t1][t2][t3];

        printf("The error between 2 random points in solution (expect 0): %f\n",value);
        printf("radiator heat at checked point %f\n",f[t1][t2][t3]);
    
    #elif TWK_OMP_GS1
        printf("hello OMP GS1!\n");
        
        
        Init_OMP2(u, u_old, f, N);
        
        ts = omp_get_wtime();
            
        gauss_seidel_OMP1(f,u,N,iter_max,&tolerance);      
            
        te = omp_get_wtime() - ts; 

        int t1 = 1;
        int t2 = 1;
        int t3 = (int) ceil((double)(N+2)/3);

        double value = -6*u[t1][t2][t3] + u[t1][t2][t3 + 1] + u[t1][t2][t3-1] + u[t1][t2-1][t3] + u[t1][t2+1][t3] + u[t1-1][t2][t3] + u[t1+1][t2][t3];
        value *= (N+1)*(N+1)/4;
        value += f[t1][t2][t3];

        printf("The error between 2 random points in solution (expect 0): %f\n",value);
        printf("radiator heat at checked point %f\n",f[t1][t2][t3]);

    #elif TWK_OMP_GS2
        printf("hello OMP GS2!\n");
        
        
        Init_OMP2(u, u_old, f, N);
        
        ts = omp_get_wtime();
            
        gauss_seidel_OMP2(f,u,N,iter_max,&tolerance);      
            
        te = omp_get_wtime() - ts;     
        
    #else
        fprintf(stderr, "MISSING ENV VARIABLE TWK_OMP{X}\n");
        return 1;
    #endif

    printf("Jacobi: Elapsed time: %f \n", te);







    // Checking Correctness
    /*

 
    


    int t1 = 1;
    int t2 = 1;
    int t3 = (int) ceil((double)(N+2)/3);

    double value = -6*u[t1][t2][t3] + u[t1][t2][t3 + 1] + u[t1][t2][t3-1] + u[t1][t2-1][t3] + u[t1][t2+1][t3] + u[t1-1][t2][t3] + u[t1+1][t2][t3];
    value *= (N+1)*(N+1)/4;
    value += f[t1][t2][t3];

    printf("The error between 2 random points in solution (expect 0): %f\n",value);
    printf("radiator heat at checked point %f\n",f[t1][t2][t3]);

    // Checking Correctness
    int t1 = 1;
    int t2 = 1;
    int t3 = (int) ceil((double)(N+2)/3);

    double value = -6*u[t1][t2][t3] + u[t1][t2][t3 + 1] + u[t1][t2][t3-1] + u[t1][t2-1][t3] + u[t1][t2+1][t3] + u[t1-1][t2][t3] + u[t1+1][t2][t3];
    value *= (N+1)*(N+1)/4;
    value += f[t1][t2][t3];

    printf("The error between 2 random points in solution (expect 0): %f\n",value);
    printf("radiator heat at checked point %f\n",f[t1][t2][t3]);
    */

    switch(output_type) {
	case 0:
	    // no output at all
	    break;
	case 3:
	    output_ext = ".bin";
	    sprintf(output_filename, "%s_%d%s", output_prefix, N, output_ext);
	    fprintf(stderr, "Write binary dump to %s: ", output_filename);
	    print_binary(output_filename, N, u);
	    break;
	case 4:
	    output_ext = ".vtk";
	    sprintf(output_filename, "%s_%d%s", output_prefix, N, output_ext);
	    fprintf(stderr, "Write VTK file to %s: ", output_filename);
	    print_vtk(output_filename, N+2, u);
	    break;
	default:
	    fprintf(stderr, "Non-supported output type!\n");
	    break;
    }

    // de-allocate memory
    free(u);
    free(f);
    free(u_old);

    return(0);
}
