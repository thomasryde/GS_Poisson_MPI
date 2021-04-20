# Makefile
#
TARGET = poisson

SOURCES	= main_OMP.c print.c alloc3d.c InitVari.c gauss_seidel.c jacobi.c
OBJECTS	= print.o alloc3d.o InitVari.o
MAIN = main_OMP.o
OBJS	= $(MAIN) gauss_seidel.o jacobi.o

# options and settings for the GCC compilers
#
CC	= gcc 
DEFS	= 
OPT	= -g -O3 -fopenmp -ffast-math #-march=broadwell
IPO	= 
ISA	= 
CHIP	= 
ARCH	= 
PARA	= 
CFLAGS	= $(DEFS) $(ARCH) $(OPT) $(ISA) $(CHIP) $(IPO) $(PARA) $(XOPTS)
LDFLAGS = -lm 

all: $(TARGET)

$(TARGET): $(OBJECTS) $(OBJS)
	$(CC) -o $@ $(CFLAGS) $(OBJS) $(OBJECTS) $(LDFLAGS)


$(MAIN):
	$(CC) -o $@ -D_GAUSS_SEIDEL -D_JACOBI $(CFLAGS) -c main_OMP.c 

clean:
	@/bin/rm -f core *.o *~

realclean: clean
	@/bin/rm -f $(TARGET)

# DO NOT DELETE

main.o: main_OMP.c print.h gauss_seidel.h jacobi.h
print.o: print.h
