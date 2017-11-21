####################################

	[Sogang Distributed programming]

	Assignment 1 - problem 6

====================================

	Usage
  ---------

  1. sequential
  
  Compile with >> mpicc seq_image_proc.c  ppm.c -o sequential.out
  run with >> mpirun -np process # -hostfile -host sequential.out ./example/Iggy.1024.ppm ./result/name.ppm

  2. parallel
  
  Compile with >> mpicc par_image_proc.c  ppm.c -o parallel.out
  run with >> mpirun -np process # -hostfile -host parallel.out ./example/Iggy.1024.ppm ./result/name.ppm


####################################
