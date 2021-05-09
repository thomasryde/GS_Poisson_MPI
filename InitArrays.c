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
        for (k = 0; k < N; k++){
            u[i][0][k] = 20;
            u[i][n-1][k] = 20;
        }
    }
    
    for (i = 0; i < n; i++){
        for (j = 0; j < n; j++){
            u[i][j][0] = 0;
            u[i][j][N-1] = 20;
        }
    }
    
    for (j = 0; j < n; j++){
        for (k = 0; k < N; k++){
            u[0][j][k] = 20;
            u[n-1][j][k] = 20;
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

void InitializeF(double ***f, int n, int N, int size, int rank){ //fucking wrong, and also should change y and z
    //double rad_x_min = -1.0;
    double rad_x_max = -3.0/8.0;
    //double rad_y_min = -1.0;
    double rad_y_max = -1.0/2.0;
    double rad_z_min = -2.0/3.0;
    double rad_z_max = 0.0;

    double eps = 2/(double)N;
    int col = rank % (int) sqrt((double) size);
    int row = rank/(int) sqrt((double) size);
    int x_start = col*(n - 2);
    int z_start = row*(n - 2);
    
    double x,y,z;
    for(int i=1; i < n-1; i++){                //x
        for(int j=1; j < n-1; j++){            //z
            for(int k=1; k < N-1; k++){        //y
                x = (x_start + i)*eps - 1;
                z = (z_start + j)*eps - 1;
                y = k*eps - 1;
                
                if (x < rad_x_max && y < rad_y_max && z < rad_z_max && z > rad_z_min){
                    f[i][j][k] = 200;
                }
                else {
                    f[i][j][k] = 0;
                }
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
    int x_start = col*(n - 2);
    int z_start = row*(n - 2);

    double x,y,z;
    double eps = 2/(double)N;

    for(int i=1; i < n-1; i++){                //x
        for(int j=1; j < n-1; j++){            //z
            for(int k=1; k < N-1; k++){        //y
                x = (x_start + i)*eps - 1;
                z = (z_start + j)*eps - 1;
                y = k*eps - 1;
                f[i][j][k] = 3*PI*PI*sin(PI*x)*sin(PI*y)*sin(PI*z);
            }
        }
    } 

  
}


void InitializeU_mp(double ***u, int n, int N){
    // Initializing u[x,z,y]
    int i,j,k;
    #pragma omp parallel
    {
    #pragma omp for
    for (i = 0; i < n; i++){
        for (k = 0; k < N; k++){
            u[i][0][k] = 20;
            u[i][n-1][k] = 20;
        }
    }
    #pragma omp for
    for (i = 0; i < n; i++){
        for (j = 0; j < n; j++){
            u[i][j][0] = 0;
            u[i][j][N-1] = 20;
        }
    }
    #pragma omp for
    for (j = 0; j < n; j++){
        for (k = 0; k < N; k++){
            u[0][j][k] = 20;
            u[n-1][j][k] = 20;
        }
    }
    #pragma omp for
    for (i = 1; i < n-1; i++){
        for (j = 1; j < n-1; j++){
            for (k = 1; k < N-1; k++){
                u[i][j][k] = 1;
            }
        }
    }
    }
}

void InitializeF_mp(double ***f, int n, int N, int size, int rank){ //fucking wrong, and also should change y and z
    //double rad_x_min = -1.0;
    double rad_x_max = -3.0/8.0;
    //double rad_y_min = -1.0;
    double rad_y_max = -1.0/2.0;
    double rad_z_min = -2.0/3.0;
    double rad_z_max = 0.0;

    double eps = 2/(double)N;
    int col = rank % (int) sqrt((double) size);
    int row = rank/(int) sqrt((double) size);
    int x_start = col*(n - 2);
    int z_start = row*(n - 2);
    
    double x,y,z;
    #pragma omp parallel for
    for(int i=1; i < n-1; i++){                //x
        for(int j=1; j < n-1; j++){            //z
            for(int k=1; k < N-1; k++){        //y
                x = (x_start + i)*eps - 1;
                z = (z_start + j)*eps - 1;
                y = k*eps - 1;
                
                if (x < rad_x_max && y < rad_y_max && z < rad_z_max && z > rad_z_min){
                    f[i][j][k] = 200;
                }
                else {
                    f[i][j][k] = 0;
                }
            }
        }
    } 
    
}