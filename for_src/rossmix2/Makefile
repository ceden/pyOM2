

include ../../site_specific.mk_${HOSTTYPE}

all: rossmix2_module.o rossmix2.o  rossmix2_parameter.o rossmix2_diag_snap.o rossmix2_friction.o

clean :
	 rm -f *.o *.mod 

rossmix2_module.o  : rossmix2_module.f90 ../main/main_module.o 
	 	$(F90) $(F90FLAGS) -I../main -c rossmix2_module.f90
rossmix2.o: rossmix2.f90 ../main/main_module.o ../etc/timing_module.o rossmix2_module.o
	 	$(F90) $(F90FLAGS)  -I../main -I../etc -I../diagnostics -c rossmix2.f90
rossmix2_parameter.o: rossmix2_parameter.f90 ../main/main_module.o ../etc/timing_module.o rossmix2_module.o
	 	$(F90) $(F90FLAGS) -I../main -I../etc  -c rossmix2_parameter.f90
rossmix2_friction.o: rossmix2_friction.f90 ../main/main_module.o ../etc/timing_module.o rossmix2_module.o
	 	$(F90) $(F90FLAGS) -I../main -I../etc  -c rossmix2_friction.f90        
rossmix2_diag_snap.o: rossmix2_diag_snap.f90 ../main/main_module.o ../etc/timing_module.o rossmix2_module.o
	 	$(F90) $(F90FLAGS) $(CDFFLAGS) -I../main -I../etc  -c rossmix2_diag_snap.f90
