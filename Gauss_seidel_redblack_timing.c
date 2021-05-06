#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "mpi.h"
#include "Checks.h"

void Gauss_seidel_redblack_timing(double ***f,double *** u, int n, int N,int max_iter,double * tolerance) {

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
    //Timing Variables
    double t1,t2;
    double time_send = 0;
    double time_compute = 0;
    

    int neigh[4];
    NeighbourCheck(neigh, size, rank);

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
    
    if (neigh[0] != -1){
        // send and recieve above
        MPI_Irecv(&(u[1][0][1]), 1, send, neigh[0], max_iter+2, MPI_COMM_WORLD,&requests[0]);
        MPI_Isend(&(u[1][1][1]), 1, send, neigh[0], max_iter+2, MPI_COMM_WORLD,&requests[1]);
    }
    if (neigh[1] != -1){
        // send and recieve left
        MPI_Irecv(&(u[0][1][0]), (n-2)*N, MPI_DOUBLE, neigh[1], max_iter+2, MPI_COMM_WORLD,&requests[2]);
        MPI_Isend(&(u[1][1][0]), (n-2)*N, MPI_DOUBLE, neigh[1], max_iter+2, MPI_COMM_WORLD,&requests[3]);
    }
    if (neigh[2] != -1){
        // send and recieve right
        MPI_Irecv(&(u[n-1][1][0]), (n-2)*N, MPI_DOUBLE, neigh[2], max_iter+2, MPI_COMM_WORLD, &requests[4]);
        MPI_Isend(&(u[n-2][1][0]), (n-2)*N, MPI_DOUBLE, neigh[2], max_iter+2, MPI_COMM_WORLD, &requests[5]);
    }
    if (neigh[3] != -1){
        // send and recieve below
        MPI_Irecv(&(u[1][n-1][1]), 1, send, neigh[3], max_iter+2, MPI_COMM_WORLD,&requests[6]);
        MPI_Isend(&(u[1][n-2][1]), 1, send, neigh[3], max_iter+2, MPI_COMM_WORLD,&requests[7]);
    }
    
    //Main Loop
    while (iter <= max_iter && FrobNorm > tolCheck){
        if (iter % FrobCheckFreq == 0){
            FrobNorm = 0;
        }
        MPI_Waitall(noRequests,requests,MPI_STATUSES_IGNORE);

        t1 = MPI_Wtime();

        for(i=1; i<n-1;i++){         //x
            for(j=1; j<n-1;j++){     //z
                for(k=1; k<N-1;k++){ //y
                    //Gauss-seidel iteration
                    
                    //Red Points
                    if( (i+j+k) % 2 == 0){
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

        time_compute += MPI_Wtime() - t1;
        t2 = MPI_Wtime();

        //Red points done, send and recv to black points
        if (neigh[0] != -1){
            // send and recieve above
            MPI_Irecv(&(u[1][0][1]), 1, send, neigh[0], iter, MPI_COMM_WORLD,&requests[0]);
            MPI_Isend(&(u[1][1][1]), 1, send, neigh[0], iter, MPI_COMM_WORLD,&requests[1]);
        }
        if (neigh[1] != -1){
            // send and recieve left
            MPI_Irecv(&(u[0][1][0]), (n-2)*N, MPI_DOUBLE, neigh[1], iter, MPI_COMM_WORLD, &requests[2]);
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

        MPI_Waitall(noRequests,requests,MPI_STATUSES_IGNORE);

        time_send += MPI_Wtime() - t2;

         for(i=1; i<n-1;i++){        //x
            for(j=1; j<n-1;j++){     //z
                for(k=1; k<N-1;k++){ //y
                    //Gauss-seidel iteration

                    //Black Points
                    if( (i+j+k) % 2 != 0) {
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

        //Black points done, send and recv for next iteration
        if (neigh[0] != -1){
            // send and recieve above
            MPI_Irecv(&(u[1][0][1]), 1, send, neigh[0], iter, MPI_COMM_WORLD, &requests[0]);
            MPI_Isend(&(u[1][1][1]), 1, send, neigh[0], iter, MPI_COMM_WORLD,&requests[1]);
        }
        if (neigh[1] != -1){
            // send and recieve left
            MPI_Irecv(&(u[0][1][0]), (n-2)*N, MPI_DOUBLE, neigh[1], iter, MPI_COMM_WORLD, &requests[2]);
            MPI_Isend(&(u[1][1][0]), (n-2)*N, MPI_DOUBLE, neigh[1], iter, MPI_COMM_WORLD,&requests[3]);
        }
        if (neigh[2] != -1){
            // send and recieve right
            MPI_Irecv(&(u[n-1][1][0]), (n-2)*N, MPI_DOUBLE, neigh[2], iter, MPI_COMM_WORLD, &requests[4]);
            MPI_Isend(&(u[n-2][1][0]), (n-2)*N, MPI_DOUBLE, neigh[2], iter, MPI_COMM_WORLD,&requests[5]);
        }
        if (neigh[3] != -1){
            // send and recieve below
            MPI_Irecv(&(u[1][n-1][1]), 1, send, neigh[3], iter, MPI_COMM_WORLD,&requests[6]);
            MPI_Isend(&(u[1][n-2][1]), 1, send, neigh[3], iter, MPI_COMM_WORLD,&requests[7]);
        }

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

    MPI_Allreduce(&time_compute,&time_compute,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&time_send,&time_send,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    if (rank == 0){
        printf("Average Compute Time = %f\n",time_compute/iter);
        printf("Average Send Time = %f\n",time_send/iter);
    }

    if (rank == 0) {
        printf("Stopped in iteration number: %d\n",iter-1);
    }
    MPI_Waitall(noRequests,requests,MPI_STATUSES_IGNORE);
    *tolerance = sqrt(FrobNorm);

}