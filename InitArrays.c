#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

void InitializeU(double ***u, int N){
    // Initializing u
    int i,j,k;
    
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
    
    for (i = 1; i < N+1; i++){
        for (j = 1; j < N+1; j++){
            for (k = 1; k < N+1; k++){
                u[i][j][k] = 1;
            }
        }
    }
}

void InitializeF(double ***f, int N){
    int i,j,k;
    double rad_i = 5/16*(N+2);
    double rad_j = 1.0/4*(N+2);
    double rad_k = 1.0/2*(N+2);
    int rad_k_init = (int) ceil((double) 1/6 *(N+2));
    
    for (i = 0; i < N+2; i++){
        for (j = 0; j < N+2; j++){
            for (k = 0; k < N+2; k++){
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
