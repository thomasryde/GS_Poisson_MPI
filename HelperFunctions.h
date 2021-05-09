

void ComputeInnerPoints(int x_interval[],int y_interval[],int z_interval[], int iter, double delta_sq,int FrobCheckFreq, double *FrobNorm, double ***f,double *** u, int redFlag);
//void ComputeInnerPoints_mp(int x_interval[],int y_interval[],int z_interval[], int iter, double delta_sq,int FrobCheckFreq, double *FrobNorm, double ***f,double *** u, int redFlag);

void SendRecieve(int neigh[],int n,int N,int iter,MPI_Request requests[],double ***u, MPI_Datatype send);
