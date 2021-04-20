#include <math.h>
#include <stdlib.h>
#include <stdio.h>

void jacobiIter(double ***f, double ***u, double ***u_old, double delta, int i, int j, int k, double h) {
    u[i][j][k] = h*(u_old[i-1][j][k] + u_old[i+1][j][k] + u_old[i][j-1][k] + u_old[i][j+1][k] + u_old[i][j][k-1] + u_old[i][j][k+1] + delta*delta*f[i][j][k]);
}
