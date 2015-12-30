#!/bin/bash                                                                                                       
#SBATCH -J ArcOn_ITER      # Job name                                                                          
#SBATCH -A Disc_Gal_Blobs #A-ph5             # Disc_Gal_Blobs  
#SBATCH -o ITER.o%j       # Name of stdout output file (%j expands to jobId)   
#SBATCH -e ITER.e%j       # Name of stdout output file (%j expands to jobId)                               
#SBATCH -p development   # Queue name                                                         
#SBATCH -n 64             #Total number of mpi tasks requested
#SBATCH -t 00:15:00      # Run time (hh:mm:ss) - 1.5 hours                                                   
#SBATCH --mail-user=michoski@ices.utexas.edu
#SBATCH --mail-type=ALL 

export MV2_ON_DEMAND_THRESHOLD=64

#ibrun -np 16 valgrind --tool=callgrind ./ArcOn -log_summary petsc_log_summary 

ibrun ./ArcOn #-log_summary petsc_log_summary -ksp_view petsc_ksp_summary
   #ibrun ./ArcOn

