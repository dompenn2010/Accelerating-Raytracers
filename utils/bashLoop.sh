for numThreads in 1 2 4 8 16 32 64 128 240
do
	OMP_NUM_THREADS=$numThreads ./rayTracing/raytracer.mic 
done
exit 0
