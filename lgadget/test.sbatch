#!/bin/bash                                                                                                                      
#SBATCH --job-name=gadget                                                                                            
#SBATCH --partition=cpuq          
#SBATCH --account=cpuq
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=yluo42@ucsc.edu 
#SBATCH --output=gadget.out                                                                                                                                                                                           
#SBATCH --time=1-00:00:00
#SBATCH --nodes=14                                                                                       
#SBATCH --ntasks-per-node=40                    
                                                                                                                                                                                      
module load openmpi 

LD_LIBRARY_PATH=/data/groups/leauthaud/yluo42/prospector/lux/gsl/lib

export LD_LIBRARY_PATH
                                                                                                                        
mpirun -np 512 ./L-Gadget2 ./parameterfiles/run.param > /data/groups/leauthaud/yluo42/prospector/lux/lgadget/parameterfiles/output_new/runcheck.txt
