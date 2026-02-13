OPTIONS=-DFFTFINE
OPTIONS+=-Dsigma8
OPTIONS+=-DPID
#OPTIONS+=-DZIPX
#OPTIONS+=-CB

MODEFILE:=$(wildcard *.f90)
OBJFILE:=$(addprefix ,$(notdir $(MODEFILE:.f90=.o)))

all: cicpower.x ic.x 
	@echo "made:" $^
ic.x: ic.o parameters.o  parameters.o pencil_fft.o# power_nu.o  variables.o
	$(FC) $(XFLAG) $(OPTIONS) $(FFTFLAG) $^ -o $@
cicpower.x: cicpower.o parameters.o pencil_fft.o powerspectrum.o# power_nu.o variables.o
	$(FC) $(XFLAG) $(OPTIONS) $(FFTFLAG) $^ -o $@

parameters.o: Makefile ../basic_functions.f08
variables.o: parameters.o
particle_initialization.o: variables.o
buffer_grid.o: variables.o
buffer_particle.o: variables.o
ic.o:parameters.o pencil_fft.o #power_nu.o
gadget_snap_fof.o:parameters.o variables.o particle_initialization.o buffer_grid.o buffer_particle.o

# cubefft.o: ../cubefft.f90
# 	$(FC) $(OFLAG) $(OPTIONS) $< -o $@ $(FFTFLAG)
# cicpower_global.o: ../neutrinos/cicpower_global.f90 variables.o cubefft.o pencil_fft.o parameters.o
# 	$(FC) $(OFLAG) $(OPTIONS) $< -o $@ $(FFTFLAG)
# power_nu.o: ../neutrinos/power_nu.f90 variables.o parameters.o cicpower_global.o
# 	$(FC) $(OFLAG) $(OPTIONS) $< -o $@ $(FFTFLAG)
parameters.o: ../parameters.f90
	$(FC) $(OFLAG) $(OPTIONS) $<
variables.o: ../variables.f90
	$(FC) $(OFLAG) $(OPTIONS) $<
particle_initialization.o: ../particle_initialization.f90
	$(FC) $(OFLAG) $(OPTIONS) $<
buffer_grid.o: ../buffer_grid.f90
	$(FC) $(OFLAG) $(OPTIONS) $<
buffer_particle.o: ../buffer_particle.f90
	$(FC) $(OFLAG) $(OPTIONS) $<
%.o: %.f90 Makefile parameters.o pencil_fft.o powerspectrum.o
	$(FC) $(OFLAG) $(OPTIONS) $(FFTFLAG) $< -o $@
pencil_fft.o: ../pencil_fft.f90 parameters.o
	$(FC) $(OFLAG) $(OPTIONS) $< $(FFTFLAG)
*.o: Makefile

clean:
	rm -vf *.mod *.o *.out *.err *.x *~
