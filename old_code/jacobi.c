/* jacobi.c - Poisson problem in 3d
 * 
 */
#include <math.h>
#include "alloc3d.h"
#include <stdlib.h>
#include <stdio.h>

void jacobi(double ***f,double ***u, double ***u_old, int N, int iter_max, double *tolerance) {
    
    //Initialize Constants
    double accum = 0.0; //accumilator for frobenius norm
    int n = 0; //jacobi iterator
    double d = 10.0;
    double delta_sq = 4.0/((N+1)*(N+1));
    double h = 1.0/6;
    int i,j,k;
    double*** tmp; // tmp pointer var for swaps

    //Main Loop
    while (n<iter_max && d > *tolerance){
        // u_old = u
        // for(i=1; i<N+1;i++){
        //     for(j=1; j<N+1;j++){
        //         for(k=1; k<N+1;k++){
        //             u_old[i][j][k] = u[i][j][k];
        //         }
        //     }
        // }

        tmp = u_old;
        u_old = u;
        u = tmp;
        
        // u = update function
        for(i=1; i<N+1;i++){
            for(j=1; j<N+1;j++){
                for(k=1; k<N+1;k++){
                    //jacobi iteration
                    u[i][j][k] = h*(u_old[i-1][j][k] + u_old[i+1][j][k] + u_old[i][j-1][k] + u_old[i][j+1][k] + u_old[i][j][k-1] + u_old[i][j][k+1] + delta_sq*f[i][j][k]);
                    
                    //computation for tolerence
                    accum += (u[i][j][k] - u_old[i][j][k]) * (u[i][j][k] - u_old[i][j][k]);
                }
            }
        }
        
        //Check Tolerence using Frobenious Norm
        d = sqrt(accum);
        accum = 0.0;
        
        n += 1;
        if(n == iter_max) {
            printf("Stopped due to iter_max \n");
        }
    }
    printf("Stopped in iteration number: %d\n",n);
    *tolerance = d;

}

void jacobi_OMP1(double ***f,double ***u, double ***u_old, int N, int iter_max, double *tolerance) {
    
    //Initialize Constants
    double accum = 0.0; //accumilator for frobenius norm
    int n = 0; //jacobi iterator
    double d = 10.0;
    double delta_sq = 4.0/((N+1)*(N+1));
    double h = 1.0/6;
    int i,j,k;
    double*** tmp; // tmp pointer var for swaps

    //Main Loop
        while (n<iter_max && d > *tolerance){
    
            tmp = u_old;
            u_old = u;
            u = tmp;
            
            #pragma omp parallel for default(none) private(i,j,k) reduction(+:accum) shared(u,u_old,f,N,delta_sq,h)
            for(i=1; i<N+1;i++){
                for(j=1; j<N+1;j++){
                    for(k=1; k<N+1;k++){
                        //jacobi iteration
                        u[i][j][k] = h*(u_old[i-1][j][k] + u_old[i+1][j][k] + u_old[i][j-1][k] + u_old[i][j+1][k] + u_old[i][j][k-1] + u_old[i][j][k+1] + delta_sq*f[i][j][k]);
                        
                        //computation for tolerence
                        accum += (u[i][j][k] - u_old[i][j][k]) * (u[i][j][k] - u_old[i][j][k]);
                    }
                }
            }
            
            
            //Check Tolerence using Frobenious Norm
            d = sqrt(accum);
            accum = 0.0;
            
            n += 1;
        //  if(n == iter_max) {
        //      printf("Stopped due to iter_max \n");
        //  }
        }
    //printf("Stopped in iteration number: %d\n",n);
    
    *tolerance = d;

}

void jacobi_OMP2(double ***f,double ***u, double ***u_old, int N, int iter_max, double *tolerance) {
    
    //Initialize Constants
    double accum = 0.0; //accumilator for frobenius norm
    int n = 0; //jacobi iterator
    double d = 10.0;
    double delta_sq = 4.0/((N+1)*(N+1));
    double h = 1.0/6;
    int i,j,k;
    double*** tmp; // tmp pointer var for swaps

    //Main Loop
    while (n<iter_max && d > *tolerance){

        tmp = u_old;
        u_old = u;
        u = tmp;
        
        #pragma omp parallel for private(i,j,k) reduction(+: accum) shared(u,u_old,f,N,delta_sq,h)
        for(i=1; i<N+1;i++){
            for(j=1; j<N+1;j++){
                for(k=1; k<N+1;k++){
                    //jacobi iteration
                    u[i][j][k] = h*(u_old[i-1][j][k] + u_old[i+1][j][k] + u_old[i][j-1][k] + u_old[i][j+1][k] + u_old[i][j][k-1] + u_old[i][j][k+1] + delta_sq*f[i][j][k]);
                    
                    //computation for tolerence
                    accum += (u[i][j][k] - u_old[i][j][k]) * (u[i][j][k] - u_old[i][j][k]);
                }
            }
        }
        
        //Check Tolerence using Frobenious Norm
        d = sqrt(accum);
        accum = 0.0;
        
        n += 1;
      //  if(n == iter_max) {
      //      printf("Stopped due to iter_max \n");
      //  }
    }
    //printf("Stopped in iteration number: %d\n",n);
    *tolerance = d;
}