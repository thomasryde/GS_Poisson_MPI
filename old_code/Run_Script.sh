EXECUTABLE=poisson

THREADS="1 2 4 8 12 16 20 24"
#THREADS="12"
for thr in $THREADS
do
    echo "---------------------------------------------------"
    echo "Running with ${thr} thread(s)!"
    OMP_NUM_THREADS=$thr OMP_PLACES=cores OMP_PROC_BIND=spread OMP_WAIT_POLICY=active ./$EXECUTABLE 100 400 0.000000001 1 0
done

