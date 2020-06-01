#!/bin/bash

#SBATCH --job-name=check_1
#SBATCH --output=check_1.out
#SBATCH --error=check_1.err
#SBATCH --partition=broadwl
#SBATCH --time=24:00:00
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=28



# --constraint=ib must be give to guarantee a job is allocated
# nodes with Infiniband
#module unload fftw2
#module unload intel
#module unload gsl
#module unload openmpi

module load intelmpi mkl intel

export I_MPI_FALLBACK=0
export I_MPI_JOB_FAST_STARTUP=enable
export I_MPI_SCALABLE_OPTIMIZATION=enable
export I_MPI_TIMER_KIND=rdtsc
#export I_MPI_DEBUG=6                                                                                                                                                  
export I_MPI_PLATFORM_CHECK=0
export I_MPI_HYDRA_PMI_CONNECT=alltoall
export I_MPI_THREAD_LEVEL_DEFAULT=FUNNELED  ##(only master thread can make MPI calls)                                                                                  

export I_MPI_FABRICS=shm:dapl
export I_MPI_DAPL_UD=off                  # you can switch UD on if you have MPI scalability troubles                                                                  
export DAPL_IB_MTU=4096
export I_MPI_PIN_DOMAIN="[FFFFFFF]"


module unload fftw2
module unload intel
module unload gsl
module unload openmpi

module load intelmpi/5.1+intel-16.0
module load fftw2/2.1.5+intelmpi-5.1+intel-16.0
module load gsl/2.2.1+intel-16.0



#module load fftw2/2.1.5+openmpi-1.6+intel-12.1
#module load gsl/2.2.1+intel-12.1


# Run the process with mpirun. Notice -n is not required. mpirun will
# automatically figure out how many processes to run from the slurm options
#mpirun -np 224 ./L-Gadget2 ./parameterfiles/test_cdm_512_150.param.1 > /scratch/midway2/susmita87/lgadget_modified_midway2/output/150_512_0.3_0.7_1_3h/runcheck.txt.4
mpirun -np 224 ./L-Gadget2 ./parameterfiles/test_cdm_512_150.param.5 > /path-to-output-directoy/runcheck.txt

#module load fftw2/2.1.5+intelmpi-5.1+intel-16.0 intelmpi/5.1+intel-16.0 gsl/2.2.1+intel-16.0
