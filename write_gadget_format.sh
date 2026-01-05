#!/bin/bash
#SBATCH -p kshcnormal
#SBATCH -N 8
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH -o write_gadget_format.out
#SBATCH -e write_gadget_format.err
module purge
module load compiler/intel/2017.5.239 mpi/intelmpi/2017.4.239
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
PROCS=$SLURM_NPROCS
EXE='./utilities/write_gadget_format.x'
export PBS_NODEFILE=`generate_pbs_nodefile`  
sort $PBS_NODEFILE | uniq -c | awk '{print $2":"$1}' > ./hostfile_$PROCS 
HOST_FILE=./hostfile_$PROCS
echo "-genvall -genv I_MPI_FABRICS=shm:ofa -machinefile $HOST_FILE -n $SLURM_NPROCS $EXE" > cafconfig.txt

$EXE
