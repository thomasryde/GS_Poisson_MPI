
default: main

MPICC = mpicc -lm -fopenmp
CFLAGS = -O3 -fopenmp -ftree-vectorize -march=native
FFLAGS = -O3 -fopenmp -ftree-vectorize -march=native

INCLUDES = -I.
OBJS = alloc3d.o
OBJS += Checks.o
OBJS += HelperFunctions.o
OBJS += InitArrays.o

alloc3d.o: alloc3d.c alloc3d.h
	$(MPICC) -c $< $(INCLUDES) $(CFLAGS)
Algorithms.o: HelperFunctions.h Gauss_Seidel_Blocked.c Gauss_Seidel_nonblocked.c Gauss_seidel_redblack.c Gauss_seidel_redblack_v3.c Gauss_seidel_redblack_v2.c Gauss_seidel_redblack_mp.c Gauss_seidel_redblack_timing.c Algorithms.h
	$(MPICC) -c $< $(INCLUDES) $(CFLAGS)
Checks.o: Checks.c Checks.h
	$(MPICC) -c $< $(INCLUDES) $(CFLAGS)
HelperFunctions.o: HelperFunctions.c HelperFunctions.h
	$(MPICC) -c $< $(INCLUDES) $(CFLAGS)
InitArrays.o: InitArrays.c InitArrays.h
	$(MPICC) -c $< $(INCLUDES) $(CFLAGS)
main.o: main.c alloc3d.h Checks.h Algorithms.h InitArrays.h HelperFunctions.h
	$(MPICC) -c $< $(INCLUDES) $(CFLAGS)
main: main.o $(OBJS)
	$(MPICC) -o $@ $^ $(LIBS)

clean:
	rm -f *.o main
