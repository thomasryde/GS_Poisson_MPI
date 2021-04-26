#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

void Gauss_Seidel_1(double ***f,double *** u, int, n, int N,int max_iter,double * tolerance){


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
    double tolCheck = tolerence*tolerence;


    MPI_Send(&(u[i][j][k]), N, MPI_DOUBLE, int dest, iter, MPI_COMM_WORLD); //Change data and destination
    MPI_Recv(&(u[i][j][k]), N, MPI_DOUBLE, int source, iter, MPI_COMM_WORLD, MPI_STATUS_IGNORE); //Change data and source
             
    int neigh[6];
    //Main Loop
    while (iter<iter_max && FrobNorm > tolCheck){
        FrobNorm = 0;
        NeighbourCheck(neigh, size, rank, iter, N);


        for(j=0; j<n;j++){         //z
            
            // send to left
            // recieve from right

            for(i=0; i<n;i++){     //x

                if(j == 0){
                    // send to up
                    // recieve from up
                } elseif(j == n-1){
                    // recieve from below
                }

                if(i == n-1){
                    // recieve from right
                }

                for(k=0; k<N;k++){ //y
                    //Gauss-seidel iteration
                    u_tmp = u[i][j][k];
                    
                    /// NEW CODE START
                    if (i == 0 || i == n-1 )

                    u[i][j][k] = h*(u[i-1][j][k] + u[i+1][j][k] + u[i][j-1][k] + u[i][j+1][k] + u[i][j][k-1] + u[i][j][k+1] + delta_sq*f[i][j][k]);
                
                    // NEW CODE END

                    //computation for tolerence
                    FrobNorm += (u[i][j][k] - u_tmp) * (u[i][j][k] - u_tmp);
                }
                if (j == n-1){
                    // send to below
                }
            }
            
            // send to right
 
        }
        
        
        
        
        iter += 1;
        if(n == iter_max) {
            printf("Stopped due to iter_max \n");
        }
    }
    printf("Stopped in iteration number: %d\n",n);
    *tolerance = sqrt(FrobNorm);

}