
include ../site_specific.mk_${HOSTTYPE}

default: with_mpi

modules:
	cd etc; make timing_module.o
	cd density; make density.o
	cd main; make main_module.o
	cd eke; make eke_module.o
	cd isoneutral; make isoneutral_module.o
	cd idemix; make idemix_module.o
	cd rossmix; make rossmix_module.o
	cd tke; make tke_module.o
	cd obc; make obc_module.o
	cd tracer; make tracer_module.o
	cd qg_filter; make qg_module.o
	cd diagnostics; make diagnostics_module.o

all:    modules
	cd main; make 
	cd advection; make
	cd tke; make
	cd eke; make
	cd idemix; make
	cd rossmix; make
	cd isoneutral; make
	cd external; make
	cd non_hydrostatic; make
	cd density; make
	cd obc; make
	cd tracer; make
	cd qg_filter; make
	cd diagnostics; make 

config.o: all config.f90 
	$(F90) $(F90FLAGS) $(CDFFLAGS) -Imain -Iisoneutral -Itke -Ieke \
	 -Iidemix -Irossmix -Idiagnostics -Idensity -Iobc -Itracer -Iqg_filter -c config.f90

with_mpi: all config.o
	cd parallel; make parallel_mpi.o
	$(F90) main/*.o advection/*.o tke/*.o eke/*.o idemix/*.o rossmix/*.o isoneutral/*.o \
	 external/*.o non_hydrostatic/*.o density/*.o obc/*.o diagnostics/*.o qg_filter/*.o \
	 tracer/*.o config.o parallel/parallel_mpi.o etc/timing_module.o\
	 $(F90FLAGS) $(MPIFLAGS) $(CDFFLAGS)-o ../bin/model.x

without_mpi: all config.o
	cd parallel; make parallel_none.o
	$(F90) main/*.o advection/*.o tke/*.o eke/*.o idemix/*.o rossmix/*.o isoneutral/*.o \
	 external/*.o non_hydrostatic/*.o density/*.o obc/*.o diagnostics/*.o qg_filter/*.o\
	 tracer/*.o config.o parallel/parallel_none.o etc/timing_module.o\
	 $(F90FLAGS) $(CDFFLAGS)-o ../bin/model.x

files = main/main_module.o main/numerics.o main/momentum_advection.o main/vertical_velocity.o 
time_step: all
	cd main; make time_step.o
	cd parallel; make parallel_mpi.o
	$(F90) $(F90FLAGS) $(CDFFLAGS) -Imain -c config.f90
	$(F90) config.o main/time_step.o parallel/parallel_mpi.o $(files) \
            -Imain  -Idensity -Iobc obc/*.o density/*.o advection/*.o $(F90FLAGS) $(MPIFLAGS)   $(CDFFLAGS)  -o ../bin/model.x
	cd main; rm time_step.o  # this is a bit dirty but otherwise two main entry points

dirs = density diagnostics eke etc external idemix rossmix isoneutral main \
       advection non_hydrostatic parallel tke obc tracer qg_filter
clean: 
	for d in $(dirs); do cd $$d && make clean && cd ..; done
	rm -f *.o *.mod 


