  ################# #XEON ######################

  # Static Scheduling  NO Recursion
  echo "---------------------------------------------------"
  echo "Xeon Static Scheduling"

  for numThreads in 1 2 4 8 16 32 64 128 256
    do
      echo $numThreads "threads"
      OMP_NUM_THREADS=$numThreads ./XeonOMP_SS_NR 
    done
  # Dynamic Scheduling NO Recursion
	
  echo "---------------------------------------------------"
  echo "Xeon Dynamic Scheduling"

  for numThreads in 1 2 4 8 16 32 64 128 256 
    do
      echo $numThreads "threads"
      OMP_NUM_THREADS=$numThreads ./XeonOMP_DS_NR
    done
  exit 0

