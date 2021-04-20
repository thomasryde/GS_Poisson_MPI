#!/bin/bash
# 02614 - High-Performance Computing, January 2018
# 
# batch script to run matmult on a dedicated server in the hpcintro
# queue
#
# Author: Bernd Dammann <bd@cc.dtu.dk>
#

#BSUB -J Amdahl_baseline_GS
#BSUB -o Amdahl_baseline_GS%J.out
#BSUB -q hpcintro
#BSUB -n 24
#BSUB -R "rusage[mem=2048]"
#BSUB -W 25

echo "Cache size info:"
echo $(lscpu | grep cache)
#lscpu

# define the driver name to use
EXECUTABLE=poisson
#./$EXECUTABLE 200 500000 1 1 0

#THREADS="1 2 4 8 12 16 20 24"
N_S="50 75 100 125 150"

for N in $N_S
do
    echo "---------------------------------------------------"
    echo "Running with N=${N}!"
    # OMP_NUM_THREADS=$thr ./$EXECUTABLE 200 5000 0.000000000001 1 0
    ./$EXECUTABLE $N 5000 0.000000000001 1 0
done

#OMP_PLACES=cores OMP_PROC_BIND=spread OMP_WAIT_POLICY=active