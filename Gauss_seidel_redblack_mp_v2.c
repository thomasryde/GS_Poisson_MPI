#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include "mpi.h"
#include "Checks.h"

void Gauss_seidel_redblack_mp_v2(double ***f,double *** u, int n, int N,int max_iter,double * tolerance) {

    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    double localNorm;

    //Initialize Constants
    double FrobNorm = 10; //accumilator for frobenius norm
    int iter = 0; //iterator
    double d = 10.0;
    double u_tmp;
    double delta_sq = 4.0/(N*N);
    double h = 1.0/6;
    int i,j,k;
    double tolCheck = (*tolerance)*(*tolerance);
    int FrobCheckFreq = 100;

    int neigh[4];
    NeighbourCheck(neigh, size, rank);

    int Two_n_2[2] = {2,n-2};
    int Two_N_2[2] = {2,N-2};
    int One_2[2] = {1,2};
    int One_n_1[2] = {1,n-1};
    int One_N_1[2] = {1,N-1};
    int n_2_n_1[2] = {n-2,n-1};
    int N_2_N_1[2] = {N-2,N-1};

    // Make stride information
    MPI_Datatype send;
    // n-2 antal blocks - N-2 antal af elementer pr block - N*n antal elementer imellem hver block
    MPI_Type_vector(n-2,N-2,N*n,MPI_DOUBLE,&send);
    MPI_Type_commit(&send);
    
    // Initialize requests as null
    int noRequests = 8;
    MPI_Request requests[noRequests];
    for(int i = 0; i < noRequests; i++){
        requests[i] = MPI_REQUEST_NULL;
    }
    

    SendRecieve(neigh,n,N,iter,requests,u,send);
    
    double t1,time;
    //Main Loop
    while (iter <= max_iter && FrobNorm > tolCheck){
        if (iter == 1) {
            MPI_Barrier(MPI_COMM_WORLD);
            if (rank == 0) {
                t1 = MPI_Wtime();
            }
        }
        if (iter % FrobCheckFreq == 0){
            FrobNorm = 0;
        }
        
        MPI_Waitall(noRequests,requests,MPI_STATUSES_IGNORE);

        //Compute RED
        #pragma omp parallel
        {
        #pragma omp master
        {
        ComputeInnerPoints(One_2,One_N_1,One_n_1,iter,delta_sq,FrobCheckFreq, &FrobNorm, f,u,0);
        ComputeInnerPoints(n_2_n_1,One_N_1,One_n_1,iter,delta_sq,FrobCheckFreq, &FrobNorm, f,u,0);
        ComputeInnerPoints(Two_n_2,One_N_1,One_2,iter,delta_sq,FrobCheckFreq, &FrobNorm, f,u,0);
        ComputeInnerPoints(Two_n_2,One_N_1,n_2_n_1,iter,delta_sq,FrobCheckFreq, &FrobNorm, f,u,0);
        ComputeInnerPoints(Two_n_2,One_2,Two_n_2,iter,delta_sq,FrobCheckFreq, &FrobNorm, f,u,0);
        ComputeInnerPoints(Two_n_2,N_2_N_1,Two_n_2,iter,delta_sq,FrobCheckFreq, &FrobNorm, f,u,0);

        SendRecieve(neigh,n,N,iter,requests,u,send);
        }
            localNorm = 0;
            #pragma omp for private(i,j,k,u_tmp) reduction(+:localNorm)
            for(i=Two_n_2[0]; i<Two_n_2[1];i++){         //x
                for(j=Two_n_2[0]; j<Two_n_2[1];j++){     //z
                    for(k=Two_N_2[0]; k<Two_N_2[1];k++){ //y
                        //Gauss-seidel iteration
                        
                        if( (i+j+k) % 2 == 0){
                            if (iter % FrobCheckFreq == 0){
                                u_tmp = u[i][j][k];
                                u[i][j][k] = h*(u[i-1][j][k] + u[i+1][j][k] + u[i][j-1][k] + u[i][j+1][k] + u[i][j][k-1] + u[i][j][k+1] + delta_sq*f[i][j][k]);
                                localNorm += (u[i][j][k] - u_tmp) * (u[i][j][k] - u_tmp);
                            } else{
                                u[i][j][k] = h*(u[i-1][j][k] + u[i+1][j][k] + u[i][j-1][k] + u[i][j+1][k] + u[i][j][k-1] + u[i][j][k+1] + delta_sq*f[i][j][k]);
                            }
                        }
                    }
                }
            }
        }
        FrobNorm += localNorm;

        MPI_Waitall(noRequests,requests,MPI_STATUSES_IGNORE);

        //Compute BLACK
        #pragma omp parallel
        {
        #pragma omp master
        {
        ComputeInnerPoints(One_2,One_N_1,One_n_1,iter,delta_sq,FrobCheckFreq, &FrobNorm, f,u,1);
        ComputeInnerPoints(n_2_n_1,One_N_1,One_n_1,iter,delta_sq,FrobCheckFreq, &FrobNorm, f,u,1);
        ComputeInnerPoints(Two_n_2,One_N_1,One_2,iter,delta_sq,FrobCheckFreq, &FrobNorm, f,u,1);
        ComputeInnerPoints(Two_n_2,One_N_1,n_2_n_1,iter,delta_sq,FrobCheckFreq, &FrobNorm, f,u,1);
        ComputeInnerPoints(Two_n_2,One_2,Two_n_2,iter,delta_sq,FrobCheckFreq, &FrobNorm, f,u,1);
        ComputeInnerPoints(Two_n_2,N_2_N_1,Two_n_2,iter,delta_sq,FrobCheckFreq, &FrobNorm, f,u,1);

        SendRecieve(neigh,n,N,iter,requests,u,send);
        }
            localNorm = 0;
            #pragma omp for private(i,j,k,u_tmp) reduction(+:localNorm)
            for(i=Two_n_2[0]; i<Two_n_2[1];i++){         //x
                for(j=Two_n_2[0]; j<Two_n_2[1];j++){     //z
                    for(k=Two_N_2[0]; k<Two_N_2[1];k++){ //y
                        //Gauss-seidel iteration
                        
                        if( (i+j+k) % 2 == 1){
                            if (iter % FrobCheckFreq == 0){
                                u_tmp = u[i][j][k];
                                u[i][j][k] = h*(u[i-1][j][k] + u[i+1][j][k] + u[i][j-1][k] + u[i][j+1][k] + u[i][j][k-1] + u[i][j][k+1] + delta_sq*f[i][j][k]);
                                FrobNorm += (u[i][j][k] - u_tmp) * (u[i][j][k] - u_tmp);
                            } else{
                                u[i][j][k] = h*(u[i-1][j][k] + u[i+1][j][k] + u[i][j-1][k] + u[i][j+1][k] + u[i][j][k-1] + u[i][j][k+1] + delta_sq*f[i][j][k]);
                            }
                        }
                    }
                }
            }

        }
        FrobNorm += localNorm;

        if (iter % FrobCheckFreq == 0) {
            MPI_Allreduce(&FrobNorm,&FrobNorm,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        }
        
        if (iter % FrobCheckFreq == 0 && rank == 0 && iter != 0){
             printf("Iter = %d --- FrobNorm = %f\n",iter,FrobNorm);
        }

        iter += 1;
        if(n == max_iter && rank == 0) {
            printf("Stopped due to iter_max \n");
        }
        
    }

    if (rank == 0) {
        printf("Stopped in iteration number: %d\n",iter-1);
    }
    MPI_Waitall(noRequests,requests,MPI_STATUSES_IGNORE);
    *tolerance = sqrt(FrobNorm);
    MPI_Type_free(&send);
    MPI_Barrier(MPI_COMM_WORLD);
    
    if (rank == 0) {
        time = MPI_Wtime() - t1;
        printf("Average iteration runtime %f seconds (skipping iter == 0) !\n",time/max_iter);
    }
}