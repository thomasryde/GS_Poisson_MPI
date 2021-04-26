
default: main

MPICC = mpicc -lm
CFLAGS = -O3 -ftree-vectorize -march=native
FFLAGS = -O3 -ftree-vectorize -march=native

INCLUDES = -I.
OBJS = alloc3d.o
OBJS += Checks.o
OBJS += Algorithms.o
OBJS += InitArrays.o

alloc3d.o: alloc3d.c alloc3d.h
	$(MPICC) -c $< $(INCLUDES) $(CFLAGS)
Checks.o: Checks.c Checks.h
	$(MPICC) -c $< $(INCLUDES) $(CFLAGS)
Algorithms.o: Algorithms.c Algorithms.h
	$(MPICC) -c $< $(INCLUDES) $(CFLAGS)
InitArrays.o: InitArrays.c InitArrays.h
	$(MPICC) -c $< $(INCLUDES) $(CFLAGS)
main.o: main.c alloc3d.h Checks.h Algorithms.h InitArrays.h
	$(MPICC) -c $< $(INCLUDES) $(CFLAGS)
main: main.o $(OBJS)
	$(MPICC) -o $@ $^ $(LIBS)

clean:
	rm -f *.o main
