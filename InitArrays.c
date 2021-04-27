#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define PI 3.14159265358979323846

void InitializeU(double ***u, int n, int N){
    // Initializing u[x,z,y]
    int i,j,k;
    
    for (i = 0; i < n; i++){
        for (k = 0; k < N; j++){
            u[i][0][k] = 20;
            u[i][n][k] = 20;
        }
    }
   
    for (i = 0; i < n; i++){
        for (j = 0; j < n; k++){
            u[i][j][0] = 0;
            u[i][j][N] = 20;
        }
    }
    
    for (j = 0; j < n; j++){
        for (k = 0; k < N; k++){
            u[0][j][k] = 20;
            u[n][j][k] = 20;
        }
    }
    
    for (i = 1; i < n-1; i++){
        for (j = 1; j < n-1; j++){
            for (k = 1; k < N-1; k++){
                u[i][j][k] = 1;
            }
        }
    }
}

void InitializeF(double ***f, int n, int N){ //fucking wrong, and also should change y and z
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
    int z_start = col*(n - 2);
    int x_start = row*(n - 2);

    double x,y,z;
    double eps = 2/(double)N;

    for(int i=1; i < n-1; i++){                //x
        for(int j=1; j < n-1; j++){            //z
            for(int k=1; k < N-1; k++){        //y
                //x = 2*(1 - eps)*(double)(x_start + i + 1)/N - 1 + eps;
                x = (x_start + i)*eps - 1;
                //z = 2*(1 - eps)*(double)(z_start + j + 1)/N - 1 + eps;
                z = (z_start + j)*eps - 1;
                //z = 1 - 
                //y = 2*(1 - eps)*(double)(k + 1)/N - 1 + eps;
                y = k*eps - 1;
                f[i][j][k] = 3*PI*PI*sin(PI*x)*sin(PI*y)*sin(PI*z);
            }
        }
    } 

  
}