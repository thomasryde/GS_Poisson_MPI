#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "mpi.h"
#include "Checks.h"

void Gauss_Seidel_nonblocked(double ***f,double *** u, int n, int N,int max_iter,double * tolerance){

    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    //Initialize Constants
    double FrobNorm = 10; //accumilator for frobenius norm
    int iter = 0; //iterator
    double u_tmp;
    double delta_sq = 4.0/(N*N);
    double h = 1.0/6;
    int i,j,k;
    double tolCheck = (*tolerance) * (*tolerance);
    int FrobCheckFreq = 100;

    int neigh[4];
    NeighbourCheck(neigh, size, rank);

    // Initialize requests as null
    int noRequests = 8;
    MPI_Request requests[noRequests];
    for(int i = 0; i < noRequests; i++){
        requests[i] = MPI_REQUEST_NULL;
    }

    //Main Loop
    while (iter <= max_iter && FrobNorm > tolCheck){
        if (iter % FrobCheckFreq == 0){
            FrobNorm = 0;
        }
        // neigh[0] = above, neigh[1] = left, neigh[2] = right, neigh[3] = below
        for (j=1; j<n-1;j++){         //z
            
            if (neigh[1] != -1){
                // send to left & recieve from left
                MPI_Isend(&u[1][j][1], N-2, MPI_DOUBLE, neigh[1], iter, MPI_COMM_WORLD,&requests[0]);
                MPI_Irecv(&u[0][j][1], N-2, MPI_DOUBLE, neigh[1], iter, MPI_COMM_WORLD,&requests[1]);
            }
            
            for (i=1; i<n-1;i++){     //x

                if (j == 1){
                    if (neigh[0] != -1){
                        // send to up
                        MPI_Isend(&u[i][1][1], N-2, MPI_DOUBLE, neigh[0], iter, MPI_COMM_WORLD,&requests[2]);
                        // recieve from up
                        MPI_Irecv(&u[i][0][1], N-2, MPI_DOUBLE, neigh[0], iter, MPI_COMM_WORLD,&requests[3]);
                    }
                    
                } else if (j == n-2){
                    if (neigh[3] != -1){
                        // recieve from below
                        MPI_Irecv(&u[i][n-1][1], N-2, MPI_DOUBLE, neigh[3], iter, MPI_COMM_WORLD,&requests[4]);
                    }
                }

                if(i == n-2){
                    if (neigh[2] != -1) {
                        // recieve from right
                        MPI_Irecv(&u[n-1][j][1], N-2, MPI_DOUBLE, neigh[2], iter, MPI_COMM_WORLD,&requests[5]);
                    }
                }

                MPI_Waitall(noRequests,requests,MPI_STATUSES_IGNORE);

                for (k=1; k<N-1;k++){ //y
                    //Gauss-seidel iteration
                    if (iter % FrobCheckFreq == 0){
                        u_tmp = u[i][j][k];
                        u[i][j][k] = h*(u[i-1][j][k] + u[i+1][j][k] + u[i][j-1][k] + u[i][j+1][k] + u[i][j][k-1] + u[i][j][k+1] + delta_sq*f[i][j][k]);
                        FrobNorm += (u[i][j][k] - u_tmp) * (u[i][j][k] - u_tmp);
                    } else{
                        u[i][j][k] = h*(u[i-1][j][k] + u[i+1][j][k] + u[i][j-1][k] + u[i][j+1][k] + u[i][j][k-1] + u[i][j][k+1] + delta_sq*f[i][j][k]);
                    }
                }
                if (j == n-2){
                    if (neigh[3] != -1) {
                        // send to below
                        MPI_Isend(&u[i][n-2][1], N-2, MPI_DOUBLE, neigh[3], iter, MPI_COMM_WORLD,&requests[6]);
                    }
                }
            }
            if (neigh[2] != -1) {
                // send to right
                MPI_Isend(&u[n-2][j][1], N-2, MPI_DOUBLE, neigh[2], iter, MPI_COMM_WORLD,&requests[7]);
            }
        }
        
        if (iter % FrobCheckFreq == 0) {
            MPI_Allreduce(&FrobNorm,&FrobNorm,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        }
        
        if (iter % FrobCheckFreq == 0 && rank == 0){
             printf("Iter = %d --- FrobNorm = %f\n",iter,FrobNorm);
        }
        iter += 1;
    }
    if (rank == 0) {
        printf("Stopped in iteration number: %d\n",iter-1);
    }
    *tolerance = sqrt(FrobNorm);
}