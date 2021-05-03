#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "mpi.h"
#include "Checks.h"

void Gauss_seidel_redblack(double ***f,double *** u, int n, int N,int max_iter,double * tolerance) {

    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

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

    // Make stride information
    MPI_Datatype send;
    // n-2 antal blocks - N-2 antal af elementer pr block - N*n antal elementer imellem hver block
    MPI_Type_vector(n-2,N-2,N*n,MPI_DOUBLE,&send);
    MPI_Type_commit(&send);
    
    if (rank == 0){
        int send_size;
        MPI_Type_size(send,&send_size);
        printf("send size %d\n", send_size/sizeof(double));
    }
    
    if (neigh[0] != -1){
        // send and recieve above
        MPI_Send(&(u[1][1][1]), 1, send, neigh[0], max_iter+2, MPI_COMM_WORLD);
        MPI_Recv(&(u[1][0][1]), 1, send, neigh[0], max_iter+2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    if (neigh[1] != -1){
        // send and recieve left
        MPI_Send(&(u[1][1][0]), (n-2)*N, MPI_DOUBLE, neigh[1], max_iter+2, MPI_COMM_WORLD);
        MPI_Recv(&(u[0][1][0]), (n-2)*N, MPI_DOUBLE, neigh[1], max_iter+2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    if (neigh[2] != -1){
        // send and recieve right
        MPI_Send(&(u[n-2][1][0]), (n-2)*N, MPI_DOUBLE, neigh[2], max_iter+2, MPI_COMM_WORLD);
        MPI_Recv(&(u[n-1][1][0]), (n-2)*N, MPI_DOUBLE, neigh[2], max_iter+2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    if (neigh[3] != -1){
        // send and recieve below
        MPI_Send(&(u[1][n-2][1]), 1, send, neigh[3], max_iter+2, MPI_COMM_WORLD);
        MPI_Recv(&(u[1][n-1][1]), 1, send, neigh[3], max_iter+2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    
    //Main Loop
    while (iter <= max_iter && FrobNorm > tolCheck){
        FrobNorm = 0;
        for(i=1; i<n-1;i++){         //x
            for(j=1; j<n-1;j++){     //z
                for(k=1; k<N-1;k++){ //y
                    //Gauss-seidel iteration
                    
                    //Red Points
                    if( (i+j+k) % 2 == 0){
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
            MPI_Send(&(u[1][1][1]), 1, send, neigh[0], iter, MPI_COMM_WORLD);
            MPI_Recv(&(u[1][0][1]), 1, send, neigh[0], iter, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        if (neigh[1] != -1){
            // send and recieve left
            MPI_Send(&(u[1][1][0]), (n-2)*N, MPI_DOUBLE, neigh[1], iter, MPI_COMM_WORLD);
            MPI_Recv(&(u[0][1][0]), (n-2)*N, MPI_DOUBLE, neigh[1], iter, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        if (neigh[2] != -1){
            // send and recieve right
            MPI_Send(&(u[n-2][1][0]), (n-2)*N, MPI_DOUBLE, neigh[2], iter, MPI_COMM_WORLD);
            MPI_Recv(&(u[n-1][1][0]), (n-2)*N, MPI_DOUBLE, neigh[2], iter, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        if (neigh[3] != -1){
            // send and recieve below
            MPI_Send(&(u[1][n-2][1]), 1, send, neigh[3], iter, MPI_COMM_WORLD);
            MPI_Recv(&(u[1][n-1][1]), 1, send, neigh[3], iter, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }


         for(i=1; i<n-1;i++){        //x
            for(j=1; j<n-1;j++){     //z
                for(k=1; k<N-1;k++){ //y
                    //Gauss-seidel iteration

                    //Black Points
                    if( (i+j+k) % 2 != 0) {
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
            MPI_Send(&(u[1][1][1]), 1, send, neigh[0], iter, MPI_COMM_WORLD);
            MPI_Recv(&(u[1][0][1]), 1, send, neigh[0], iter, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        if (neigh[1] != -1){
            // send and recieve left
            MPI_Send(&(u[1][1][0]), (n-2)*N, MPI_DOUBLE, neigh[1], iter, MPI_COMM_WORLD);
            MPI_Recv(&(u[0][1][0]), (n-2)*N, MPI_DOUBLE, neigh[1], iter, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        if (neigh[2] != -1){
            // send and recieve right
            MPI_Send(&(u[n-2][1][0]), (n-2)*N, MPI_DOUBLE, neigh[2], iter, MPI_COMM_WORLD);
            MPI_Recv(&(u[n-1][1][0]), (n-2)*N, MPI_DOUBLE, neigh[2], iter, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        if (neigh[3] != -1){
            // send and recieve below
            MPI_Send(&(u[1][n-2][1]), 1, send, neigh[3], iter, MPI_COMM_WORLD);
            MPI_Recv(&(u[1][n-1][1]), 1, send, neigh[3], iter, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        if (iter % FrobCheckFreq == 0) {
            MPI_Allreduce(&FrobNorm,&FrobNorm,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        }
        
        if (iter % FrobCheckFreq == 0 && rank == 0){
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
    *tolerance = sqrt(FrobNorm);

}