for numThreads in 1 2 4 8 16 32 64 128 256
   do
           OMP_NUM_THREADS=$numThreads ./raytracer
   done
   exit 0
