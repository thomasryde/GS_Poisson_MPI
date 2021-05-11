// main.c - Poisson problem in 3D

#include <stdio.h>
#include <stdlib.h>
#include "alloc3d.h"
#include <math.h>
#include "InitArrays.h"
#include <omp.h>
#include <mpi.h>
#include "HelperFunctions.h"
#include "Algorithms.h"
#include "Checks.h"


#define N_DEFAULT 100

//#define case_type 4 // 0 for correctness test,
                      // 1 for blocked Gauss Seidel,
                      // 2 for non blocked Gauss Seidel, 
                      // 3 for Red and Black Gauss Seidel.
                      // 4 for Red and Black Gauss Seidel && OpenMP.
                      // 5 for Red and Black Gauss Seidel v2.
                      // 6 for Red and Black Gauss Seidel v3.
                      // 7 for Red and Black Gauss Seidel && OpenMP NUMA
                      // 8 for Red and Black Gauss Seidel && OpenMPv2.
                      // 100 for Red and Black Gauss Seidel with Timing.

int main(int argc, char *argv[]) {

    MPI_Init(&argc, &argv);
    
    int size, rank;

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int 	N = N_DEFAULT;
    int 	iter_max = 1000;
    double 	tolerance;
    int	    case_type;
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
    case_type   = atof(argv[4]);  // start T for all inner grid points
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
        //exit(-1);
        free(u);

        MPI_Finalize();
        return 0;
    }
    if ( (f = d_malloc_3d(n, n, N)) == NULL ) {
        printf("array f: allocation failed on rank %d\n",rank);
        //exit(-1);
        free(u);
        free(f);

        MPI_Finalize();
        return 0;
    }

    switch(case_type){
    case 0: //Correctness test
        if (rank == 0){
            printf("---Running correctness check---\n");
        }
        InitU_cortest(u, n, N);
        InitF_cortest(f, n, N, size, rank);

        MPI_Barrier(MPI_COMM_WORLD);
        Gauss_Seidel_Blocked(f,u,n,N,iter_max,&tolerance);
        if ( (u_anal = d_malloc_3d(n, n, N)) == NULL ) {
            printf("array u_anal: allocation failed on rank %d\n",rank);
            exit(-1);
        }

        double error = CorrectnessCheck(u,u_anal,rank,int_sqrt_size,n,N);
        free(u_anal);
        break;

    case 1: // Blocked Gauss Seidel
        if (rank == 0){
            printf("---Running blocked Gauss Seidel---\n");
        }
        
        InitializeU(u, n, N);
        InitializeF(f, n, N, size, rank);
        
        Gauss_Seidel_Blocked(f,u,n,N,iter_max,&tolerance);
        
        break;

    case 2: // non-blocked Gauss Seidel
        if (rank == 0){
            printf("---Running non blocked Gauss Seidel---\n");
        }
        InitializeU(u, n, N);
        InitializeF(f, n, N, size, rank);
        
        Gauss_Seidel_nonblocked(f,u,n,N,iter_max,&tolerance);
        
        break;

    case 3: // Red and Black Gauss Seidel non blocked
        if (rank == 0){
            printf("---Running Gauss Seidel with red and black---\n");
        }
        InitializeU(u, n, N);
        InitializeF(f, n, N, size, rank);
        
        Gauss_seidel_redblack(f,u,n,N,iter_max,&tolerance);
        
        break;
    case 4: // Red and Black Gauss Seidel with openMP
        if (rank == 0){
            printf("---Running HYBRID Gauss Seidel with red and black and OpenMP---\n");
        }
        InitializeU(u, n, N);
        InitializeF(f, n, N, size, rank);
        
        Gauss_seidel_redblack_mp(f,u,n,N,iter_max,&tolerance);
        
        break;
    case 5: // Red and Black Gauss Seidel non blocked v2!!
        if (rank == 0){
            printf("---Running Gauss Seidel with red and black v2---\n");
        }
        InitializeU(u, n, N);
        InitializeF(f, n, N, size, rank);
        
        Gauss_seidel_redblack_v2(f,u,n,N,iter_max,&tolerance);
        
        break;
    case 6: // Red and Black Gauss Seidel non blocked v3!!
        if (rank == 0){
            printf("---Running Gauss Seidel with red and black v3---\n");
        }
        InitializeU(u, n, N);
        InitializeF(f, n, N, size, rank);
        
        Gauss_seidel_redblack_v3(f,u,n,N,iter_max,&tolerance);
        
        break;
    case 7: // Red and Black Gauss Seidel with openMP NUMA
        if (rank == 0){
            printf("---Running HYBRID Gauss Seidel with red and black and OpenMP NUMA---\n");
        }
        InitializeU_mp(u, n, N);
        InitializeF_mp(f, n, N, size, rank);
        
        Gauss_seidel_redblack_mp(f,u,n,N,iter_max,&tolerance);
        
        break;
    case 8: // Red and Black Gauss Seidel with openMP_v2
        if (rank == 0){
            printf("---Running HYBRID Gauss Seidel with red and black and OpenMP V2---\n");
        }
        InitializeU_mp(u, n, N);
        InitializeF_mp(f, n, N, size, rank);
        
        Gauss_seidel_redblack_mp_v2(f,u,n,N,iter_max,&tolerance);
        
        break;
    case 100: // Red and Black Gauss Seidel non blocked with timing
        if (rank == 0){
            printf("---Running Gauss Seidel TIMING with red and black---\n");
        }
        InitializeU(u, n, N);
        InitializeF(f, n, N, size, rank);
        MPI_Barrier(MPI_COMM_WORLD);
        if (rank == 0){
            t1 = MPI_Wtime();
        }
        Gauss_seidel_redblack_timing(f,u,n,N,iter_max,&tolerance);
        MPI_Barrier(MPI_COMM_WORLD);
        if (rank == 0){
            time = MPI_Wtime() - t1;
            printf("Average iteration runtime %f seconds!\n",time/iter_max);
        }
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
    for (int i = 0; i<4; i++) {
        printf("Rank %d neigbour[%d] %d\n",rank,i,neigh[i]);
    }
    */

    // de-allocate memory
    free(u);
    free(f);

    MPI_Finalize();

    return 0;
    
}
