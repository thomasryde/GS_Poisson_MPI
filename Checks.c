#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PI 3.14159265358979323846

int InputCheck(int N, int size){
    int int_sqrt_size;
    if (sqrt((double)size) * sqrt((double)size) == size){
        int_sqrt_size = (int) sqrt((double) size);
    } else{
        fprintf(stderr, "You have entered an invalid size!");
        exit(-1);
    }

    if (N % int_sqrt_size != 0){
        fprintf(stderr, "Grid size does not match with number of processors! Exiting...\n");
        exit(-1);
    }
    return int_sqrt_size;
}


double CorrectnessCheck(double *** u, double *** u_anal, int rank,int int_sqrt_size,int n,int N){
    
    if (rank == 0){
        printf("Make sure to use the correct f and u boundary conditions!!\n");
    }

    double FrobError = 0;
    
    int col = rank % int_sqrt_size;
    int row = rank/int_sqrt_size;
    

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
                u_anal[i][j][k] = sin(PI*x)*sin(PI*y)*sin(PI*z);
                FrobError += (u_anal[i][j][k] - u[i][j][k])*(u_anal[i][j][k] - u[i][j][k]);
            }
        }
    }        

    return FrobError;

}




int CheckEdge(int size, int rank) {
    //edge above = 1, edge left = 2, edge right = 3, edge below = 4,
    int edge, len;

    len = (int) sqrt((double) size);
    edge = -1;
    
    if (rank < len) {
        edge = 1;
    }
    if ((rank % len) == 0){
        edge = 2;
    }
    if (((rank+1) % len) == 0){
        edge = 3;
    }
    if (rank + 1 > size-len) {
        edge = 4;
    }
    if ((rank < len) && ((rank % len) == 0)) {
        edge = 12;
    }
    if ((rank < len) && (((rank+1) % len) == 0)) {
        edge = 13;
    }
    if (((rank % len) == 0) && (rank + 1 > size-len)) {
        edge = 24;
    }
    if ((((rank+1) % len) == 0) && (rank + 1 > size-len)) {
        edge = 34;
    }

    return edge;
}

void NeighbourCheck(int neigh[], int size, int rank, int iter, int N){
    // neigh[0] = above, neigh[1] = left, neigh[2] = right, neigh[3] = below, neigh[4] = in front?, neigh[5] = behind?
    int edge1, edge2,len;
    len = (int) sqrt(size);
    
    neigh[0] = rank - len;
    neigh[1] = rank - 1;
    neigh[2] = rank + 1;
    neigh[3] = rank + len;
    neigh[4] = rank;
    neigh[5] = rank;
    
    int edge = CheckEdge(size, rank);
    
    if (edge != -1){
        neigh[edge-1] = -1;
    }
    if (iter == 0){
        neigh[4] = -1;
    }
    if (iter == (N-1)){
        neigh[5] = -1;
    }
    if (edge > 10) {
        edge1 = edge/10;
        edge2 = edge % 10;
        neigh[edge1 - 1] = -1;
        neigh[edge2 - 1] = -1;
    }
}

