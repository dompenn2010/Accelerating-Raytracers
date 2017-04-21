# OpenACC

  ################## KERNELS ####################
  # No Recursion No RegCountMax
  
  echo "---------------------------------------------"
  echo "OpenACC NoCountMax Kernels"
    
  for numSpheres in 32 1000 5000 10000 20000 30000 40000
    do
      echo $numSpheres "spheres" 
      ./Kernels_NR_NRC $numSpheres
    done

  ################## Collapsed Loops ###################

  # No Recursion No RegCountMax
  echo "---------------------------------------------"
  echo "OpenACC NoCountMax Collapsed Loops"

  for numSpheres in 32 1000 5000 10000 20000 30000 40000
    do
      echo $numSpheres "spheres" 
      ./Collapsed_NR_NRC $numSpheres
    done
  exit 0

