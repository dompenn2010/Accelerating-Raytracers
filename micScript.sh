  ################## MIC #######################
  # Static Scheduling  No Recursion
  echo "------------------------------------------------"
  echo "MIC Static Scheduling"
	
  for numThreads in 50 75 100 125 150 200 250
    do
      echo $numThreads "threads" 
      OMP_NUM_THREADS=$numThreads ./MicOMP_SS_NR.mic
    done

  echo "------------------------------------------------"
  echo "MIC Dynamic Scheduling"

  for numThreads in 50 75 100 125 150 200 250
    do
      echo $numThreads "threads" 
      OMP_NUM_THREADS=$numThreads ./MicOMP_DS_NR.mic
    done
  exit 0

