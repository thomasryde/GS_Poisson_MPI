#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"

void ComputeInnerPoints(int x_interval[],int y_interval[],int z_interval[], int iter, double delta_sq,int FrobCheckFreq, double *FrobNorm, double ***f,double *** u, int redFlag){
    double u_tmp;
    double h = 1.0/6;
    for(int i=x_interval[0]; i<x_interval[1];i++){         //x
        for(int j=z_interval[0]; j<z_interval[1];j++){     //z
            for(int k=y_interval[0]; k<y_interval[1];k++){ //y
                //Gauss-seidel iteration
                
                //Red Points
                if( (i+j+k) % 2 == redFlag){
                    if (iter % FrobCheckFreq == 0){
                        u_tmp = u[i][j][k];
                        u[i][j][k] = h*(u[i-1][j][k] + u[i+1][j][k] + u[i][j-1][k] + u[i][j+1][k] + u[i][j][k-1] + u[i][j][k+1] + delta_sq*f[i][j][k]);
                        *FrobNorm += (u[i][j][k] - u_tmp) * (u[i][j][k] - u_tmp);
                    } else{
                        u[i][j][k] = h*(u[i-1][j][k] + u[i+1][j][k] + u[i][j-1][k] + u[i][j+1][k] + u[i][j][k-1] + u[i][j][k+1] + delta_sq*f[i][j][k]);
                    }
                }
            }
        }
    }
}

void SendRecieve(int neigh[],int n,int N,int iter,MPI_Request requests[],double ***u, MPI_Datatype send){
    if (neigh[0] != -1){
        // send and recieve above
        MPI_Irecv(&(u[1][0][1]), 1, send, neigh[0], iter, MPI_COMM_WORLD,&requests[0]);
        MPI_Isend(&(u[1][1][1]), 1, send, neigh[0], iter, MPI_COMM_WORLD,&requests[1]);
    }
    if (neigh[1] != -1){
        // send and recieve left
        MPI_Irecv(&(u[0][1][0]), (n-2)*N, MPI_DOUBLE, neigh[1], iter, MPI_COMM_WORLD,&requests[2]);
        MPI_Isend(&(u[1][1][0]), (n-2)*N, MPI_DOUBLE, neigh[1], iter, MPI_COMM_WORLD,&requests[3]);
    }
    if (neigh[2] != -1){
        // send and recieve right
        MPI_Irecv(&(u[n-1][1][0]), (n-2)*N, MPI_DOUBLE, neigh[2], iter, MPI_COMM_WORLD,&requests[4]); 
        MPI_Isend(&(u[n-2][1][0]), (n-2)*N, MPI_DOUBLE, neigh[2], iter, MPI_COMM_WORLD,&requests[5]);
    }
    if (neigh[3] != -1){
        // send and recieve below
        MPI_Irecv(&(u[1][n-1][1]), 1, send, neigh[3], iter, MPI_COMM_WORLD,&requests[6]);
        MPI_Isend(&(u[1][n-2][1]), 1, send, neigh[3], iter, MPI_COMM_WORLD,&requests[7]);
    }
}

// void ComputeInnerPoints_mp(int x_interval[],int y_interval[],int z_interval[], int iter, double delta_sq,int FrobCheckFreq, double *FrobNorm, double ***f,double *** u, int redFlag){
//     double u_tmp;
//     double h = 1.0/6;
//     int i,j,k;
//     static double localNorm = 0;
//     #pragma omp for private(i,j,k,u_tmp) reduction(+:localNorm)
//     for(i=x_interval[0]; i<x_interval[1];i++){         //x
//         for(j=z_interval[0]; j<z_interval[1];j++){     //z
//             for(k=y_interval[0]; k<y_interval[1];k++){ //y
//                 //Gauss-seidel iteration
                
//                 //Red Points
//                 if( (i+j+k) % 2 == redFlag){
//                     if (iter % FrobCheckFreq == 0){
//                         u_tmp = u[i][j][k];
//                         u[i][j][k] = h*(u[i-1][j][k] + u[i+1][j][k] + u[i][j-1][k] + u[i][j+1][k] + u[i][j][k-1] + u[i][j][k+1] + delta_sq*f[i][j][k]);
//                         FrobNorm += (u[i][j][k] - u_tmp) * (u[i][j][k] - u_tmp);
//                     } else{
//                         u[i][j][k] = h*(u[i-1][j][k] + u[i+1][j][k] + u[i][j-1][k] + u[i][j+1][k] + u[i][j][k-1] + u[i][j][k+1] + delta_sq*f[i][j][k]);
//                     }
//                 }
//             }
//         }
//     }
// }