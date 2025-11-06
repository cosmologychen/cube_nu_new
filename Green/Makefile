MODEFILE:=$(wildcard *.f90)
OBJFILE:=$(addprefix ,$(notdir $(MODEFILE:.f90=.o)))
#OPTIONS+=-Wall

all: kernel_initialization.x force_test_pairwise.x force_test_NFW.x
	@echo "made:" $^
kernel_initialization.x: kernel_initialization.o parameters.o Green.o
	$(FC) $(XFLAG) $(OPTIONS) $^ -o $@
force_test_pairwise.x: force_test_pairwise.o parameters.o Green.o
	$(FC) $(XFLAG) $(OPTIONS) $(FFTFLAG) $^ -o $@
force_test_NFW.x: force_test_NFW.o parameters.o variables.o Green.o kick.o apm.o pp.o buffer_grid.o buffer_particle.o
	$(FC) $(XFLAG) $(OPTIONS) $(FFTFLAG) $^ -o $@

force_test_NFW.o: variables.o pp.o
kick.o: variables.o pp.o

parameters.o: ../parameters.f90
	$(FC) $(OFLAG) $(OPTIONS) $<
buffer_grid.o: ../main/buffer_grid.f90 parameters.o variables.o
	$(FC) $(OFLAG) $(OPTIONS) $< -o $@
buffer_particle.o: ../main/buffer_particle.f90 parameters.o variables.o
	$(FC) $(OFLAG) $(OPTIONS) $< -o $@

%.o: %.f90 Makefile parameters.o variables.o #pencil_fft.o cubefft.o
	$(FC) $(OFLAG) $(OPTIONS) $(FFTFLAG) $< -o $@

*.o: Makefile

clean:
	rm -f *.mod *.o *.out *.x *~
