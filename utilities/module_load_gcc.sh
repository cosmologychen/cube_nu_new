module purge
module load compiler/devtoolset/7.3.1 mpi/hpcx/2.7.4/gcc-7.3.1 mathlib/fftw/3.3.8/single/gnu
export PATH=$PATH:/public/home/acd4q6s2ve/software/OpenCoarrays-hpcx/bin/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/public/home/acd4q6s2ve/software/OpenCoarrays-hpcx/lib64/

export FC='caf'
export XFLAG_NO_OMP='-O3 -cpp -mcmodel=large'
export XFLAG='-O3 -cpp -fopenmp -mcmodel=large'
export OFLAG_NO_OMP=$XFLAG_NO_OMP' -c'
export OFLAG=$XFLAG' -c'
export FFTDIR='/public/software/mathlib/fftw/3.3.8/single/gnu/'
export FFTFLAG='-I'$FFTDIR'include/ -L'$FFTDIR'lib/ -lfftw3f -lfftw3f_threads -lm -ldl'

export OMP_STACKSIZE=16000M
#export KMP_STACKSIZE=16000M
export OMP_NUM_THREADS=32
#export OMP_THREAD_LIMIT=4
export FOR_COARRAY_NUM_IMAGES=8
ulimit
ulimit -s unlimited

# run executable by ./a.out and set number of images by FOR_COARRAY_NUM_IMAGES
