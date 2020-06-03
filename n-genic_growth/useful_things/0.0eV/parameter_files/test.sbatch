#!/bin/bash                                                                                                                      
#SBATCH --job-name=test_CDM                                                                                            
#SBATCH --partition=defq           
#SBATCH --mail-type=FAIL
##SBATCH --mail-user=yluo42@ucsc.edu 
#SBATCH --output=test_CDM.out                                                                                             
#SBATCH --error=test_CDM.err                                                                                              
#SBATCH --time=00:45:00
#SBATCH --nodes=8                                                                                                       
#SBATCH --ntasks-per-node=16                                                                                            
#SBATCH --output=test_CDM.log
                                                                                                                                                                                       
module load openmpi                                                                                                                                                                                                             
LD_LIBRARY_PATH=/data/groups/leauthaud/yluo42/prospector/lux/gsl/lib
export LD_LIBRARY_PATH
                                                                                                                        
mpirun -np 128 /data/groups/leauthaud/yluo42/prospector/lux/n-genic_growth/N-GenIC N-GenIC.param
