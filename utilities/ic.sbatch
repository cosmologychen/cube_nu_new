#!/bin/bash
#SBATCH -J initial_conditions
#SBATCH -p kshctest
#SBATCH -N 1
#SBATCH --ntasks-per-node=1

#export LD_LIBRARY_PATH=/opt/gridview/slurm17/lib/:$LD_LIBRARY_PATH
#export I_MPI_PMI_LIBRARY=/opt/gridview/slurm17/lib/libpmi.so
#export I_MPI_PMI_LIBRARY=/opt/gridview/slurm17/lib/libpmi2.so


export OMP_STACKSIZE=16000M
#export KMP_STACKSIZE=16000M
export OMP_NUM_THREADS=1
#export OMP_THREAD_LIMIT=4
export FOR_COARRAY_NUM_IMAGES=1
ulimit
ulimit -s unlimited

srun ./ic.x
