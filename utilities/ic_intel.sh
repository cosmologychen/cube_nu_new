#!/bin/bash
#SBATCH -p kshctest
#SBATCH -N 8
#SBATCH --ntasks-per-node=1
#SBATCH -o %j.out
#SBATCH -e %j.err
module purge
module load compiler/intel/2017.5.239 mpi/intelmpi/2017.4.239
export PATH=/public/home/acd4q6s2ve/test:$PATH
PROCS=$SLURM_NPROCS
srun hostname|sort|uniq -c|awk '{print $2":"$1}' > hostsfile
#srun -n 8 hostname > hostsfile
#export OMP_NUM_THREADS=8
#export FOR_COARRAY_NUM_IMAGES=8
./ic.x
