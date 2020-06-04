#!/bin/bash                                                                                                                      
#SBATCH --job-name=LCDM                                                                                            
#SBATCH --partition=defq           
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=yluo42@ucsc.edu 
#SBATCH --output=LCDM.out                                                                                             
#SBATCH --error=LCDM.err                                                                                              
#SBATCH --time=72:00:00
#SBATCH --nodes=32                                                                                         
#SBATCH --ntasks-per-node=16                    
#SBATCH --output=LCDM.log
                                                                                                                                                                                      
module load openmpi 

LD_LIBRARY_PATH=/data/groups/leauthaud/yluo42/prospector/lux/gsl/lib

export LD_LIBRARY_PATH

                                                                                                                        
mpirun -np 512 ./L-Gadget2 ./parameterfiles/run.param > /data/groups/leauthaud/yluo42/prospector/lux/lgadget/parameterfiles/output/runcheck.txt
