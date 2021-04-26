#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define PI 3.14159265358979323846

void InitializeU(double ***u, int n, int N){
    // Initializing u
    int i,j,k;
    /*
    for (i = 0; i < N+2; i++){
        for (j = 0; j < N+2; j++){
            u[i][j][0] = 20;
            u[i][j][N+1] = 20;
        }
    }
   
    for (i = 0; i < N+2; i++){
        for (k = 0; k < N+2; k++){
            u[i][0][k] = 0;
            u[i][N+1][k] = 20;
        }
    }
    
    for (j = 0; j < N+2; j++){
        for (k = 0; k < N+2; k++){
            u[0][j][k] = 20;
            u[N+1][j][k] = 20;
        }
    }
    */
    for (i = 0; i < n; i++){
        for (j = 0; j < n; j++){
            for (k = 0; k < N; k++){
                u[i][j][k] = 1;
            }
        }
    }
}

void InitializeF(double ***f, int n, int N){
    int i,j,k;
    double rad_i = 5/16*(n);
    double rad_j = 1.0/4*(n);
    double rad_k = 1.0/2*(N);
    int rad_k_init = (int) ceil((double) 1/6 *(N));
    
    for (i = 0; i < n; i++){
        for (j = 0; j < n; j++){
            for (k = 0; k < N; k++){
                f[i][j][k] = 0;
            }
        }
    }

    for (i = 0; i <= rad_i; i++){
        for (j = 0; j <= rad_j; j++){
            for (k = rad_k_init; k <= rad_k; k++){
                    f[i][j][k] = 200;
            }
        }
    }
}


void InitU_cortest(double ***u, int n, int N){
    int i,j,k;
  
    for (i = 0; i < n; i++){
        for (j = 0; j < n; j++){
            for (k = 0; k < N; k++){
                u[i][j][k] = 0;
            }
        }
    }
}

void InitF_cortest(double ***f, int n, int N, int size, int rank) {
    int i,j,k;
    
    int col = rank % (int) sqrt((double) size);
    int row = rank/(int) sqrt((double) size);
    int z_start = col*n;
    int x_start = row*n;

    double x,y,z;
    double eps = (double)2/N;

    for(int i=0; i < n; i++){                //x
        for(int j=0; j < n; j++){            //z
            for(int k=0; k < N; k++){        //y
                x = 2*(1 - eps)*(double)(x_start + i + 1)/N - 1 + eps;
                z = 2*(1 - eps)*(double)(z_start + j + 1)/N - 1 + eps;
                y = 2*(1 - eps)*(double)(k + 1)/N - 1 + eps;
                f[i][j][k] = 3*PI*PI*sin(PI*x)*sin(PI*y)*sin(PI*z);
            }
        }
    } 

  
}