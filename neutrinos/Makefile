OPTIONS=-DFFTFINE
OPTIONS+=-Dsigma8
OPTIONS+=-DPID
XFLAG=-O3 -fpp -qopenmp  -coarray=distributed  -mcmodel=large
OFLAG=-O3 -fpp -qopenmp  -coarray=distributed  -mcmodel=large -c
#OPTIONS+=-CB

MODEFILE:=$(wildcard *.f90)
OBJFILE:=$(addprefix ,$(notdir $(MODEFILE:.f90=.o)))

all: test.x
	@echo "made:" $^
test.x: variables.o cubefft.o pencil_fft.o powerspectrum.o parameters.o power_nu.o test.o cicpower_segment.o
	$(FC) $(XFLAG) $(OPTIONS) $(FFTFLAG) $^ -o $@

parameters.o: Makefile ../basic_functions.f08
variables.o: parameters.o
power_nu.o: variables.o parameters.o cicpower_segment.o
test.o: power_nu.o variables.o parameters.o power_nu.o

parameters.o: ../parameters.f90
	$(FC) $(OFLAG) $(OPTIONS) $<
variables.o: ../variables.f90
	$(FC) $(OFLAG) $(OPTIONS) $<
pencil_fft.o: ../pencil_fft.f90 parameters.o
	$(FC) $(OFLAG) $(OPTIONS) $< $(FFTFLAG)
cubefft.o: ../cubefft.f90 parameters.o variables.o
	$(FC) $(OFLAG) $(OPTIONS) $< -o $@ $(FFTFLAG)
powerspectrum.o: ../utilities/powerspectrum.f90  pencil_fft.o
	$(FC) $(OFLAG) $(OPTIONS) $< -o $@ $(FFTFLAG)
cicpower_segment.o: cicpower_segment.f90 cubefft.o pencil_fft.o powerspectrum.o parameters.o variables.o
	$(FC) $(OFLAG) $(OPTIONS) $< -o $@ $(FFTFLAG)
power_nu.o: power_nu.f90 variables.o parameters.o cicpower_segment.o
	$(FC) $(OFLAG) $(OPTIONS) $<
test.o: test.f90 power_nu.o variables.o parameters.o
	$(FC) $(OFLAG) $(OPTIONS) $<
*.o: Makefile

clean:
	rm -vf *.mod *.o *.out *.err *.x *~
