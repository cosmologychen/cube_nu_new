module purge
module load compiler/intel/2017.5.239 mpi/intelmpi/2017.4.239 mathlib/fftw/3.3.8/single/intel

export FC='ifort'
export XFLAG_NO_OMP='-O3 -fpp -coarray=distributed -mcmodel=large -coarray-config-file=cafconfig.txt'
export XFLAG='-O3 -fpp -qopenmp -coarray=distributed -mcmodel=large -coarray-config-file=cafconfig.txt'
export OFLAG_NO_OMP=$XFLAG_NO_OMP' -c'
export OFLAG=$XFLAG' -c'
#export FFTFLAG='-I/public/software/compiler/intel/intel-compiler-2017.5.239/mkl/include/fftw/ -mkl'
export FFTDIR='/public/software/mathlib/fftw/3.3.8/single/intel/'
export FFTFLAG='-I'$FFTDIR'include/ -L'$FFTDIR'lib/ -lfftw3f -lfftw3f_threads -lm -ldl'

export OMP_STACKSIZE=16000M
#export KMP_STACKSIZE=16000M
export OMP_NUM_THREADS=32
#export OMP_THREAD_LIMIT=4
export FOR_COARRAY_NUM_IMAGES=8
ulimit
ulimit -s unlimited

# run executable by ./a.out and set number of images by FOR_COARRAY_NUM_IMAGES
