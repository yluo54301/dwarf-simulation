#!/bin/bash                                                                                                                                                                                                                  
#SBATCH --partition=cpuq           
#SBATCH --account=cpuq
#SBATCH --job-name=LCDM
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=yluo42@ucsc.edu 
#SBATCH --output=LCDM.out                                                                                                                                                                                           
#SBATCH --time=03:00:00
#SBATCH --ntasks=16                                                                                              
#SBATCH --nodes=1                                                                                         
#SBATCH --ntasks-per-node=16                    
                                                                                                                                                                                      
module load openmpi 

LD_LIBRARY_PATH=/data/groups/leauthaud/yluo42/prospector/lux/gsl/lib

export LD_LIBRARY_PATH

                                                                                                                        
mpirun -np 16 ./L-Gadget2 ./parameterfiles/run.param > /data/groups/leauthaud/yluo42/prospector/lux/lgadget/runcheck.txt
