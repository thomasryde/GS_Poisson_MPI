#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "mpi.h"
#include "Checks.h"

void Gauss_Seidel_Blocked(double ***f,double *** u, int n, int N,int max_iter,double * tolerance){
    
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

    // Neighbors
    int neigh[4];
    NeighbourCheck(neigh, size, rank);
    
    double t1,time;
    //Main Loop
    while (iter <= max_iter && FrobNorm > tolCheck){
        //----------- TIMING
        if (iter == 1) {
            MPI_Barrier(MPI_COMM_WORLD);
            if (rank == 0) {
                t1 = MPI_Wtime();
            }
        }
        if (iter % FrobCheckFreq == 0){
            FrobNorm = 0;
        }

        // neigh[0] = above, neigh[1] = left, neigh[2] = right, neigh[3] = below
        for (j=1; j<n-1;j++){         //z
            
            if (neigh[1] != -1){
                // send to left & recieve from left
                MPI_Send(&u[1][j][1], N-2, MPI_DOUBLE, neigh[1], iter, MPI_COMM_WORLD);
                MPI_Recv(&u[0][j][1], N-2, MPI_DOUBLE, neigh[1], iter, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            
            for (i=1; i<n-1;i++){     //x

                if (j == 1){
                    if (neigh[0] != -1){
                        // send to up
                        MPI_Send(&u[i][1][1], N-2, MPI_DOUBLE, neigh[0], iter, MPI_COMM_WORLD);
                        // recieve from up
                        MPI_Recv(&u[i][0][1], N-2, MPI_DOUBLE, neigh[0], iter, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    }
                    
                } else if (j == n-2){
                    if (neigh[3] != -1){
                        // recieve from below
                        MPI_Recv(&u[i][n-1][1], N-2, MPI_DOUBLE, neigh[3], iter, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    }
                }

                if(i == n-2){
                    if (neigh[2] != -1) {
                        // recieve from right
                        MPI_Recv(&u[n-1][j][1], N-2, MPI_DOUBLE, neigh[2], iter, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    }
                }
                
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
                        MPI_Send(&u[i][n-2][1], N-2, MPI_DOUBLE, neigh[3], iter, MPI_COMM_WORLD);
                    }
                }
            }
            if (neigh[2] != -1) {
                // send to right
                MPI_Send(&u[n-2][j][1], N-2, MPI_DOUBLE, neigh[2], iter, MPI_COMM_WORLD);
            }
        }
        
        // Compute FrobNorm
        if (iter % FrobCheckFreq == 0) {
            MPI_Allreduce(&FrobNorm,&FrobNorm,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        }
        
        if (iter % FrobCheckFreq == 0 && rank == 0 && iter != 0){
             printf("Iter = %d --- FrobNorm = %f\n",iter,FrobNorm);
        }
        iter += 1;
    }
    if (rank == 0) {
        printf("Stopped in iteration number: %d\n",iter-1);
    }
    *tolerance = sqrt(FrobNorm);
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) {
        time = MPI_Wtime() - t1;
        printf("Average iteration runtime %f seconds (skipping iter == 0) !\n",time/max_iter);
    }
}