/* gauss_seidel.c - Poisson problem in 3d
 *
 */
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

void gauss_seidel(double ***f,double ***u, int N, int iter_max, double *tolerance) {
    //Initialize Constants
    double accum = 0.0; //accumilator for frobenius norm
    int n = 0; //jacobi iterator
    double d = 10.0;
    double u_tmp;
    double delta_sq = 4.0/((N+1)*(N+1));
    double h = 1.0/6;
    int i,j,k;

    //Main Loop
    while (n<iter_max && d > *tolerance){
        for(i=1; i<N+1;i++){
            for(j=1; j<N+1;j++){
                for(k=1; k<N+1;k++){
                    //Gauss-seidel iteration
                    u_tmp = u[i][j][k];
                    u[i][j][k] = h*(u[i-1][j][k] + u[i+1][j][k] + u[i][j-1][k] + u[i][j+1][k] + u[i][j][k-1] + u[i][j][k+1] + delta_sq*f[i][j][k]);
                    
                    //computation for tolerence
                    accum += (u[i][j][k] - u_tmp) * (u[i][j][k] - u_tmp);
                }
            }
        }
        
        //Check Tolerence using Frobenious Norm
        d = sqrt(accum);
        accum = 0.0;
        /*
        if(n%100 == 0) {
            printf("%f \n",d);
        }
        */
        
        n += 1;
        if(n == iter_max) {
            printf("Stopped due to iter_max \n");
        }
    }
    printf("Stopped in iteration number: %d\n",n);
    *tolerance = d;
}


void gauss_seidel_OMP_base(double ***f,double ***u, int N, int iter_max, double *tolerance) {
    //Initialize Constants
    double accum = 0.0; //accumilator for frobenius norm
    int n = 0; //jacobi iterator
    double d = 10.0;
    double u_tmp;
    double delta_sq = 4.0/((N+1)*(N+1));
    double h = 1.0/6;
    int i,j,k;

    //Main Loop
    while (n<iter_max && d > *tolerance){
        #pragma omp parallel for private(i,j,k,u_tmp) reduction(+: accum) shared(u,f,N,delta_sq,h)
        for(i=1; i<N+1;i++){
            for(j=1; j<N+1;j++){
                for(k=1; k<N+1;k++){
                    //Gauss-seidel iteration
                    u_tmp = u[i][j][k];
                    u[i][j][k] = h*(u[i-1][j][k] + u[i+1][j][k] + u[i][j-1][k] + u[i][j+1][k] + u[i][j][k-1] + u[i][j][k+1] + delta_sq*f[i][j][k]);
                    //computation for tolerence
                    accum += (u[i][j][k] - u_tmp) * (u[i][j][k] - u_tmp);
                }
            }
        }
        
        //Check Tolerence using Frobenious Norm
        d = sqrt(accum);
        accum = 0.0;
        
        /*
        if(n%100 == 0) {
            printf("%f \n",d);
        }
        */
       
        
        n += 1;
        if(n == iter_max) {
            printf("Stopped due to iter_max \n");
        }
    }
    printf("Stopped in iteration number: %d\n",n);
    *tolerance = d;
}

void gauss_seidel_OMP1(double ***f,double ***u, int N, int iter_max, double *tolerance) {
    //Initialize Constants
    double accum = 0.0; //accumilator for frobenius norm
    int n = 0; //jacobi iterator
    double d = 10.0;
    double u_tmp;
    double delta_sq = 4.0/((N+1)*(N+1));
    double h = 1.0/6;
    int i,j,k;

    //Main Loop
    while (n<iter_max && d > *tolerance){
        #pragma omp parallel private(i,j,k,u_tmp) reduction(+: accum) shared(u,f,N,delta_sq,h) 
        {
        #pragma omp for ordered(3) schedule(static,1)
        for(i=1; i<N+1;i++){
            for(j=1; j<N+1;j++){
                for(k=1; k<N+1;k++){
                    #pragma omp ordered depend(sink: i-1,j,k) depend(sink: i,j-1,k) depend(sink: i,j,k-1)
                    //Gauss-seidel iteration
                    u_tmp = u[i][j][k];
                    u[i][j][k] = h*(u[i-1][j][k] + u[i+1][j][k] + u[i][j-1][k] + u[i][j+1][k] + u[i][j][k-1] + u[i][j][k+1] + delta_sq*f[i][j][k]);
                    //computation for tolerence
                    accum += (u[i][j][k] - u_tmp) * (u[i][j][k] - u_tmp);
                    #pragma omp ordered depend(source)
                } 
            }
        }
        } // end of parallelization
        
        //Check Tolerence using Frobenious Norm
        d = sqrt(accum);
        accum = 0.0;
        
        /*
        if(n%100 == 0) {
            printf("%f \n",d);
        }
        */
        
        n += 1;
        if(n == iter_max) {
            printf("Stopped due to iter_max \n");
        }
    }
    printf("Stopped in iteration number: %d\n",n);
    *tolerance = d;
}

void gauss_seidel_OMP2(double ***f,double ***u, int N, int iter_max, double *tolerance) {
    //Initialize Constants
    double accum = 0.0; //accumilator for frobenius norm
    int n = 0; //jacobi iterator
    double d = 10.0;
    double u_tmp;
    double delta_sq = 4.0/((N+1)*(N+1));
    double h = 1.0/6;
    int i,j,k;

    //Main Loop
    while (n<iter_max && d > *tolerance){
        #pragma omp parallel private(i,j,k,u_tmp) reduction(+: accum) shared(u,f,N,delta_sq,h) 
        {
        #pragma omp for ordered(2) schedule(static,1)
        for(i=1; i<N+1;i++){
            for(j=1; j<N+1;j++){
            #pragma omp ordered depend(sink: i-1,j) depend(sink: i,j-1)
                for(k=1; k<N+1;k++){

                    //Gauss-seidel iteration
                    u_tmp = u[i][j][k];
                    u[i][j][k] = h*(u[i-1][j][k] + u[i+1][j][k] + u[i][j-1][k] + u[i][j+1][k] + u[i][j][k-1] + u[i][j][k+1] + delta_sq*f[i][j][k]);
                    //computation for tolerence

                    accum += (u[i][j][k] - u_tmp) * (u[i][j][k] - u_tmp);
                } 
            #pragma omp ordered depend(source)
            }
        }
        } // end of parallelization
        
        //Check Tolerence using Frobenious Norm
        d = sqrt(accum);
        accum = 0.0;
        
        /*
        if(n%100 == 0) {
            printf("%f \n",d);
        }
        */
        
        n += 1;
        if(n == iter_max) {
            printf("Stopped due to iter_max \n");
        }
    }
    printf("Stopped in iteration number: %d\n",n);
    *tolerance = d;
}
