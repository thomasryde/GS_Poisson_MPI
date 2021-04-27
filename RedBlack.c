#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "mpi.h"
#include "Checks.h"

void Gauss_seidel_redblack(double ***f,double *** u, int, n, int N,int max_iter,double * tolerance) {

    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    //Initialize Constants
    double FrobNorm = 0; //accumilator for frobenius norm
    int iter = 0; //iterator
    double d = 10.0;
    double u_tmp;
    double delta_sq = 4.0/((N+1)*(N+1));
    double h = 1.0/6;
    int i,j,k;
    double tolCheck = (*tolerence)*(*tolerence);
       
    int neigh[6];
    NeighbourCheck(neigh, size, rank);

    //Main Loop
    while (iter<max_iter && FrobNorm > tolCheck){
        FrobNorm = 0;
        for(i=1; i<n-1;i++){         //x
            for(j=1; j<n-1;j++){     //z
                for(k=1; k<N-1;k++){ //y
                    //Gauss-seidel iteration
                    
                    //Red Points
                    if( (j+k) % 2 == 0){
                        u_tmp = u[i][j][k];

                        u[i][j][k] = h*(u[i-1][j][k] + u[i+1][j][k] + u[i][j-1][k] + u[i][j+1][k] + u[i][j][k-1] + u[i][j][k+1] + delta_sq*f[i][j][k]);
                        
                        //computation for tolerence
                        FrobNorm += (u[i][j][k] - u_tmp) * (u[i][j][k] - u_tmp);
                    }
                    
                }
            }
        }

        //Red points done, send and recv to black points
        if (neigh[0] != -1){
            // send and recieve above
            MPI_Send(&(u[1][n-1][1]), (n-2)*(N-1), MPI_DOUBLE, neigh[0], iter, MPI_COMM_WORLD);
            MPI_Recv(&(u[0][n-1][1]), (n-2)*(N-1), MPI_DOUBLE, neigh[0], iter, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        if (neigh[1] != -1){
            // send and recieve left
            MPI_Send(&(u[1][1][1]), (n-2)*(N-1), MPI_DOUBLE, neigh[1], iter, MPI_COMM_WORLD);
            MPI_Recv(&(u[0][0][1]), (n-2)*(N-1), MPI_DOUBLE, neigh[1], iter, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        if (neigh[2] != -1){
            // send and recieve right
            MPI_Send(&(u[n-1][1][1]), (n-2)*(N-1), MPI_DOUBLE, neigh[2], iter, MPI_COMM_WORLD);
            MPI_Recv(&(u[n][0][1]), (n-2)*(N-1), MPI_DOUBLE, neigh[2], iter, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        if (neigh[3] != -1){
            // send and recieve below
            MPI_Send(&(u[1][n-1][1]), (n-2)*(N-1), MPI_DOUBLE, neigh[3], iter, MPI_COMM_WORLD);
            MPI_Recv(&(u[0][n][1]), (n-2)*(N-1), MPI_DOUBLE, neigh[3], iter, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }


        for(i=1; i<N-1;i++){         //y
            for(j=1; j<n-1;j++){     //z
                for(k=1; k<n-1;k++){ //x
                    //Gauss-seidel iteration

                    //Black Points
                    if( (j+k) % 2 != 0) {
                        u_tmp = u[i][j][k];

                        u[i][j][k] = h*(u[i-1][j][k] + u[i+1][j][k] + u[i][j-1][k] + u[i][j+1][k] + u[i][j][k-1] + u[i][j][k+1] + delta_sq*f[i][j][k]);
                        
                        //computation for tolerence
                        FrobNorm += (u[i][j][k] - u_tmp) * (u[i][j][k] - u_tmp);
                    }
                    
                    
                }
            }
        }
        //Black points done, send and recv for next iteration
        if (neigh[0] != -1){
            // send and recieve above
            MPI_Send(&(u[1][j][1]), n-2, MPI_DOUBLE, neigh[0], iter, MPI_COMM_WORLD);
            MPI_Recv(&(u[0][j][1]), n-2, MPI_DOUBLE, neigh[0], iter, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        if (neigh[1] != -1){
            // send and recieve left
            MPI_Send(&(u[1][j][1]), n-2, MPI_DOUBLE, neigh[1], iter, MPI_COMM_WORLD);
            MPI_Recv(&(u[0][j][1]), n-2, MPI_DOUBLE, neigh[1], iter, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        if (neigh[2] != -1){
            // send and recieve right
            MPI_Send(&(u[1][j][1]), n-2, MPI_DOUBLE, neigh[2], iter, MPI_COMM_WORLD);
            MPI_Recv(&(u[0][j][1]), n-2, MPI_DOUBLE, neigh[2], iter, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        if (neigh[3] != -1){
            // send and recieve below
            MPI_Send(&(u[1][j][1]), n-2, MPI_DOUBLE, neigh[3], iter, MPI_COMM_WORLD);
            MPI_Recv(&(u[0][j][1]), n-2, MPI_DOUBLE, neigh[3], iter, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        iter += 1;
        if(n == iter_max) {
            printf("Stopped due to iter_max \n");
        }
    }
    printf("Stopped in iteration number: %d\n",n);
    *tolerance = sqrt(FrobNorm);

}