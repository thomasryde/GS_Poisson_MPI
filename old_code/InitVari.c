
#include <math.h>
#include "alloc3d.h"
#include <stdlib.h>
#include <stdio.h>
#include <omp.h>

void InitU(double ***u, int N){
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

void InitF(double ***f, int N){
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


void Init_OMP2(double ***u, double ***u_old, double ***f, int N){
    // Initializing u
    int i,j,k;
    double rad_j = 1/4*(N+2);
    double rad_k = 1/2*(N+2);
    int rad_k_init = (int) ceil((double) 1/6 *(N+2));

    #pragma omp parallel default(none) private(i,j,k) shared(u,u_old,f,N,rad_j,rad_k,rad_k_init)
    {
    #pragma omp for nowait
    for (i = 0; i < N+2; i++){
        for (j = 0; j < N+2; j++){
            u[i][j][0] = 20;
            u[i][j][N+1] = 20;
        }
    }
    
    #pragma omp for nowait
    for (i = 0; i < N+2; i++){
        for (k = 0; k < N+2; k++){
            u[i][0][k] = 0;
            u[i][N+1][k] = 20;
        }
    }
    
    #pragma omp for nowait
    for (j = 0; j < N+2; j++){
        for (k = 0; k < N+2; k++){
            u[0][j][k] = 20;
            u[N+1][j][k] = 20;
        }
    }
    #pragma omp for nowait
    for (i = 1; i < N+1; i++){
        for (j = 1; j < N+1; j++){
            for (k = 1; k < N+1; k++){
                u[i][j][k] = 1;
            }
        }
    }

    // INIT U_OLD
    #pragma omp for nowait
    for (i = 0; i < N+2; i++){
        for (j = 0; j < N+2; j++){
            u_old[i][j][0] = 20;
            u_old[i][j][N+1] = 20;
        }
    }
    #pragma omp for nowait
    for (i = 0; i < N+2; i++){
        for (k = 0; k < N+2; k++){
            u_old[i][0][k] = 0;
            u_old[i][N+1][k] = 20;
        }
    }
    #pragma omp for nowait
    for (j = 0; j < N+2; j++){
        for (k = 0; k < N+2; k++){
            u_old[0][j][k] = 20;
            u_old[N+1][j][k] = 20;
        }
    }
    #pragma omp for nowait
    for (i = 1; i < N+1; i++){
        for (j = 1; j < N+1; j++){
            for (k = 1; k < N+1; k++){
                u_old[i][j][k] = 1;
            }
        }
    }

    // // Init F
    // #pragma omp for nowait
    // for (i = 0; i < N+2; i++){
    //     for (j = 0; j < N+2; j++){
    //         for (k = 0; k < N+2; k++){
    //             f[i][j][k] = 0;
                
    //         }
    //     }
    // }

    // INIT F
    #pragma omp for nowait
    for (i = 0; i < N+2; i++){
        for (j = 0; j < N+2; j++){
            for (k = 0; k < N+2; k++){
                f[i][j][k] = 0;
                
            }
        }
    }
    
    #pragma omp for nowait
    for (i = 0; i < 5/16*(N+2); i++){
         for (j = 0; j <= rad_j; j++){
             for (k = rad_k_init; k <= rad_k; k++){
                f[i][j][k] = 200;
             }
         }
    }

    } // end of pragma
}
