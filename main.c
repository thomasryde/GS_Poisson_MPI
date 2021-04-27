// main.c - Poisson problem in 3D

#include <stdio.h>
#include <stdlib.h>
#include "alloc3d.h"
#include <math.h>
#include "InitArrays.h"
#include <omp.h>

#include <mpi.h>
#include "Algorithms.h"
#include "Checks.h"

#define N_DEFAULT 100

#define case_type 1 // 0 for timing, 1 for correctness test.

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
    double 	***u_anal = NULL;
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

    int int_sqrt_size = InputCheck(N,size);

    int n = N/int_sqrt_size + 2;
    N = N + 2;
    double t1,time;

    // allocate memory
    if ( (u = d_malloc_3d(n, n, N)) == NULL ) {
        printf("array u: allocation failed on rank %d\n",rank);
        exit(-1);
    }
    if ( (f = d_malloc_3d(n, n, N)) == NULL ) {
        printf("array f: allocation failed on rank %d\n",rank);
        exit(-1);
    }
    
    switch(case_type){
    case 0: // Timing case
        InitializeU(u, n, N);
        InitializeF(f, n, N);

        MPI_Barrier(MPI_COMM_WORLD);
        if (rank == 0){
            t1 = MPI_Wtime();
        }
        MPI_Barrier(MPI_COMM_WORLD);
        if (rank == 0){
            time = MPI_Wtime() - t1;
            printf("It took %f seconds!\n",time);
        }
        break;

    case 1: //Correctness test
        InitU_cortest(u, n, N);
        InitF_cortest(f, n, N,size,rank);

        Gauss_Seidel_1(f,u,n,N,iter_max,&tolerance);
        
        if ( (u_anal = d_malloc_3d(n, n, N)) == NULL ) {
            printf("array u_anal: allocation failed on rank %d\n",rank);
            exit(-1);
        }
        double error = CorrectnessCheck(u,u_anal,rank,int_sqrt_size,n,N);
        printf("Error: %f, rank: %d\n",error,rank);
        break;
    default:
        printf("Incorrect case_type, chose 0 for timing - 1 for correctness check\n");
        break;
    }
    /*
    switch(output_type) {
	case 0:
	    // no output at all
	    break;
	case 3:
	    output_ext = ".bin";
	    sprintf(output_filename, "%s_%d%s", output_prefix, N, output_ext);
	    fprintf(stderr, "Write binary dump to %s: ", output_filename);
	    print_binary(output_filename, n, u);
	    break;
	case 4:
	    output_ext = ".vtk";
	    sprintf(output_filename, "%s_%d%s", output_prefix, N, output_ext);
	    fprintf(stderr, "Write VTK file to %s: ", output_filename);
	    print_vtk(output_filename, n, u);
	    break;
	default:
	    fprintf(stderr, "Non-supported output type!\n");
	    break;
    }
    */

    /* // this works :smiley:
    printf("Does rank %d have a edge %d\n",rank,CheckEdge(size,rank)); 
    int neigh[6];
    NeighbourCheck(neigh, size, rank, 10, N);
    for (int i = 0; i<6; i++) {
        printf("Rank %d neigbour[%d] %d\n",rank,i,neigh[i]);
    }
    */

    // de-allocate memory
    free(u);
    free(f);

    MPI_Finalize();

    return 0;
    
}
