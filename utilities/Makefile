OPTIONS=-DFFTFINE
OPTIONS+=-Dsigma8
OPTIONS+=-DPID
OPTIONS+=-DZIPX
#OPTIONS+=-CB

MODEFILE:=$(wildcard *.f90)
OBJFILE:=$(addprefix ,$(notdir $(MODEFILE:.f90=.o)))

all: cicpower.x ic.x gadget_snap_fof.x fof.x dsp.x # fof_merge.x #write_gadget_format.x
	@echo "made:" $^
ic.x: ic.o parameters.o  parameters.o pencil_fft.o# power_nu.o  variables.o
	$(FC) $(XFLAG) $(OPTIONS) $(FFTFLAG) $^ -o $@
cicpower.x: cicpower.o parameters.o pencil_fft.o powerspectrum.o# power_nu.o variables.o
	$(FC) $(XFLAG) $(OPTIONS) $(FFTFLAG) $^ -o $@
dsp.x: parameters.o pencil_fft.o powerspectrum.o displacement.o
	$(FC) $(XFLAG) $(OPTIONS) $(FFTFLAG) $^ -o $@
write_gadget_format.x: write_gadget_format.o parameters.o
	$(FC) $(XFLAG) $(OPTIONS) $^ -o $@
fof.x: parameters.o variables.o particle_initialization.o buffer_grid.o buffer_particle.o fof.o
	$(FC) $(XFLAG) $(OPTIONS) $^ -o $@
fof_merge.x: parameters.o variables.o fof_merge.o
	$(FC) $(XFLAG) $(OPTIONS) $^ -o $@
gadget_snap_fof.x: parameters.o variables.o gadget_snap_fof.o particle_initialization.o buffer_grid.o buffer_particle.o
	$(FC) $(XFLAG) $(OPTIONS) $^ -o $@
multifields.x: parameters.o multifields.o
	$(FC) $(XFLAG) $(OPTIONS) $^ -o $@

parameters.o: Makefile ../basic_functions.f08
variables.o: parameters.o
particle_initialization.o: variables.o
buffer_grid.o: variables.o
buffer_particle.o: variables.o
ic.o:parameters.o pencil_fft.o #power_nu.o
gadget_snap_fof.o:parameters.o variables.o particle_initialization.o buffer_grid.o buffer_particle.o

# cubefft.o: ../cubefft.f90
# 	$(FC) $(OFLAG) $(OPTIONS) $< -o $@ $(FFTFLAG)
# cicpower_segment.o: ../neutrinos/cicpower_segment.f90 variables.o cubefft.o pencil_fft.o parameters.o
# 	$(FC) $(OFLAG) $(OPTIONS) $< -o $@ $(FFTFLAG)
# power_nu.o: ../neutrinos/power_nu.f90 variables.o parameters.o cicpower_segment.o
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
