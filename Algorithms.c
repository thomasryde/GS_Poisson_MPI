#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "mpi.h"
#include "Checks.h"

void Gauss_Seidel_1(double ***f,double *** u, int n, int N,int max_iter,double * tolerance){

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

    int neigh[4];
    NeighbourCheck(neigh, size, rank);
    
    //Main Loop
    while (iter < max_iter && FrobNorm > tolCheck){
        FrobNorm = 0;
        // neigh[0] = above, neigh[1] = left, neigh[2] = right, neigh[3] = below
        MPI_Barrier(MPI_COMM_WORLD);
        
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
                    u_tmp = u[i][j][k];
                    u[i][j][k] = h*(u[i-1][j][k] + u[i+1][j][k] + u[i][j-1][k] + u[i][j+1][k] + u[i][j][k-1] + u[i][j][k+1] + delta_sq*f[i][j][k]);
                    //computation for tolerence
                    //if(iter%500==0 && i==3 && j == 25 && k == 25){
                    //    printf("iter: %d \t rank: %d\n u_tmp = %f\n u[25][25][25] = %f\n",iter, rank, u_tmp, u[i][j][k]);
                    //}
                    FrobNorm += (u[i][j][k] - u_tmp) * (u[i][j][k] - u_tmp);
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
        
        MPI_Allreduce(&FrobNorm,&FrobNorm,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        iter += 1;
        if (iter % 100 == 0 && rank == 0){
             printf("FrobNorm = %f\n",FrobNorm);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        
    }
    printf("Stopped in iteration number: %d\n",iter);
    *tolerance = sqrt(FrobNorm);
    if (rank == 0) {
        printf("u[25][25][25] = %f\n",u[25][25][25]);
    }
}