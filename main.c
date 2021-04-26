// main.c - Poisson problem in 3D

#include <stdio.h>
#include <stdlib.h>
#include "alloc3d.h"
#include "print.h"
#include <math.h>
#include "InitVari.h"
#include <omp.h>

#include <mpi.h>
#include "gauss_seidel.h"

#define N_DEFAULT 100

int main(int argc, char *argv[]) {

    MPI_Init(&argc, &argv);
    
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int 	N = N_DEFAULT;
    int 	iter_max = 1000;
    double 	tolerance;
    double	start_T;
    int		output_type = 0;
    char	*output_prefix = "poisson_res";
    char    *output_ext    = "";
    char	output_filename[FILENAME_MAX];
    double 	***u = NULL;
    double 	***f = NULL;

    /* get the paramters from the command line */
    N         = atoi(argv[1]);	// grid size
    iter_max  = atoi(argv[2]);  // max. no. of iterations
    tolerance = atof(argv[3]);  // tolerance
    start_T   = atof(argv[4]);  // start T for all inner grid points
    if (argc == 6) {
	    output_type = atoi(argv[5]);  // output type
    }

    InputCheck(N,size);

    // allocate memory
    if ( (u = d_malloc_3d(N+2, N+2, N+2)) == NULL ) {
        perror("array u: allocation failed on rank %d\n",rank);
        exit(-1);
    }
    if ( (f = d_malloc_3d(N+2, N+2, N+2)) == NULL ) {
        perror("array f: allocation failed on rank %d\n",rank;
        exit(-1);
    }
    printf("luder");
    InitializeU(u, N);
    InitializeF(f, N);

    MPI_BARRIER();
    if (rank == 0){
        double t1 = MPI_Wtime();
    }
    Gauss_Seidel_1(f,u,N,iter_max,&tolerance);
    MPI_BARRIER();
    if (rank == 0){
        double time = MPI_Wtime() - t1;
        printf("It took %f seconds!\n",time);
    }


    if ( (u_anal = d_malloc_3d(N+2, N+2, N+2)) == NULL ) {
        perror("array u_anal: allocation failed on rank %d\n",rank);
        exit(-1);
    }
    CorrectnessCheck(double *** u, double *** u_anal,int rank, int size);

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

    MPI_Finalize();

    return(0);
    
}
