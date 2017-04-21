# Comment Frenzy

# Serial Code Here

	./serialRaytracer
# OMP

  ################## XEON ######################

  # Static Scheduling  NO Recursion
	
  for numThreads in 1 2 4 8 16 32 64 128 256
    do 
      OMP_NUM_THREADS=$numthreads ./XeonOMP_SS_NR
    done
  exit 0

  # Dynamic Scheduling NO Recursion

  for numThreads in 1 2 4 8 16 32 64 128 256
    do 
      OMP_NUM_THREADS=$numthreads ./XeonOMP_DS_NR
    done
  exit 0

  # Static Scheduling WITH Recursion
	
  for numThreads in 1 2 4 8 16 32 64 128 256
    do 
      OMP_NUM_THREADS=$numthreads ./XeonOMP_SS_WR
    done
  exit 0

  # Dynamic Scheduling WITH Recursion

  for numThreads in 1 2 4 8 16 32 64 128 256
    do 
      OMP_NUM_THREADS=$numthreads ./XeonOMP_DS_WR
    done
  exit 0


  ################## MIC #######################
  # Static Scheduling  No Recursion
	
  for numThreads in 1 2 4 8 16 32 64 128 256
    do 
      OMP_NUM_THREADS=$numthreads ./MicOMP_SS_NR
    done

  # Dynamic Scheduling No Recursion

  for numThreads in 1 2 4 8 16 32 64 128 256
    do 
      OMP_NUM_THREADS=$numthreads ./MicOMP_DS_NR
    done
  
#  # Static Scheduling WITH Recursion
#	
#  for numThreads in 1 2 4 8 16 32 64 128 256
#    do 
#      OMP_NUM_THREADS=$numthreads ./MicOMP_SS_WR
#    done
#
#  # Dynamic Scheduling WITH Recursion
#
#  for numThreads in 1 2 4 8 16 32 64 128 256
#    do 
#      OMP_NUM_THREADS=$numthreads ./MicOMP_DS_WR
#    done
  exit 0

# OpenACC

  ################## KERNELS ####################
  # No Recursion No RegCountMax

  for numSpheres in 32 1000 5000 10000 20000 30000 40000
    do 
      ./Kernels_NR_NRC $numSpheres
    done
  exit 0

  # No Recursion WITH RegCountMax

  for numSpheres in 32 1000 5000 10000 20000 30000 40000
    do 
      ./Kernels_NR_WRC $numSpheres
    done
  exit 0

  # WITH Recursion No RegCountMax

  for numSpheres in 32 1000 5000 10000 20000 30000 40000
    do 
      ./Kernels_WR_NRC $numSpheres
    done
  exit 0

  # WITH Recursion WITH RegCountMax

  for numSpheres in 32 1000 5000 10000 20000 30000 40000
    do 
      ./Kernels_WR_WRC $numSpheres
    done
  exit 0

  ################## Collapsed Loops ###################

  # No Recursion No RegCountMax

  for numSpheres in 32 1000 5000 10000 20000 30000 40000
    do 
      ./Collapsed_NR_NRC $numSpheres
    done
  exit 0

  # No Recursion WITH RegCountMax

  for numSpheres in 32 1000 5000 10000 20000 30000 40000
    do 
      ./Collapsed_NR_WRC $numSpheres
    done
  exit 0

  # WITH Recursion No RegCountMax

  for numSpheres in 32 1000 5000 10000 20000 30000 40000
    do 
      ./Collapsed_WR_NRC $numSpheres
    done
  exit 0

  # WITH Recursion WITH RegCountMax

  for numSpheres in 32 1000 5000 10000 20000 30000 40000
    do 
      ./Collapsed_WR_WRC $numSpheres
    done
  exit 0

