#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void InputCheck(int N, int size){
    if (sqrt(size) * sqrt(size) != size){
     fprintf(stderr, "You have entered an invalid size!");
     exit(-1);   
    }
    
    if (N % sqrt(size) != 0){
        fprintf(stderr, "Grid size does not match with number of processors! Exiting...\n");
        exit(-1);
    }
}


void CorrectnessCheck(u,u_anal,rank,size){
    printf("Make sure to use the correct f and u boundary conditions!\n");
    // u = sin(pi*x)*sin(pi*y)*sin(pi*z)

}