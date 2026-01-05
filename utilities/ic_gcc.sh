#!/bin/bash
#SBATCH -N 8
#SBATCH --ntasks-per-node=1
#SBATCH -p kshctest
#SBATCH --mem 0 

module purge
module load compiler/devtoolset/7.3.1 mpi/hpcx/2.7.4/gcc-7.3.1 mathlib/fftw/3.3.8/single/gnu
export FOR_COARRAY_NUM_IMAGES=$SLURM_NTASKS
export PATH=$PATH:/public/home/acd4q6s2ve/software/OpenCoarrays-hpcx/bin/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/public/home/acd4q6s2ve/software/OpenCoarrays-hpcx/lib64/

export I_MPI_FABRICS=shm:dapl
export I_MPI_DAPL_UD=enable
export I_MPI_FALLBACK_DEVICE=disable
export I_MPI_DEBUG=0
export I_MPI_PIN=disable
export I_MPI_ADJUST_REDUCE=2
export I_MPI_ADJUST_ALLREDUCE=2
export I_MPI_ADJUST_BCAST=1
export I_MPI_PLATFORM=auto
export I_MPI_DAPL_SCALABLE_PROGRESS=1
export I_MPI_DAPL_UD_PROVIDER=ofa-v2-mlx5_0-1u
export I_MPI_PIN_DOMAIN=core

cafrun -n $SLURM_NTASKS ./ic.x
