
include ../../site_specific.mk_${HOSTTYPE}

all: rossmix_module.o rossmix.o rossmix_integrate.o rossmix_reflect.o rossmix_group_velocity.o \
     rossmix_flux.o rossmix_friction.o 

clean :
	 rm -f *.o *.mod 

rossmix_module.o  : rossmix_module.f90 ../main/main_module.o 
	 	$(F90) $(F90FLAGS) -I../main -c rossmix_module.f90
rossmix.o: rossmix.f90 ../main/main_module.o ../etc/timing_module.o rossmix_module.o
	 	$(F90) $(F90FLAGS) -I../main -I../etc  -c rossmix.f90
rossmix_flux.o: rossmix_flux.f90 ../main/main_module.o rossmix_module.o
	 	$(F90) $(F90FLAGS) -I../main -c rossmix_flux.f90
rossmix_integrate.o: rossmix_integrate.f90 ../main/main_module.o ../etc/timing_module.o rossmix_module.o
	 	$(F90) $(F90FLAGS) -I../main -I../etc -c rossmix_integrate.f90
rossmix_friction.o: rossmix_friction.f90 ../main/main_module.o rossmix_module.o
	 	$(F90) $(F90FLAGS) -I../main -c rossmix_friction.f90
rossmix_reflect.o: rossmix_reflect.f90 ../main/main_module.o rossmix_module.o
	 	$(F90) $(F90FLAGS) -I../main -c rossmix_reflect.f90
rossmix_group_velocity.o: rossmix_group_velocity.f90 ../main/main_module.o rossmix_module.o
	 	$(F90) $(F90FLAGS) -I../main -c rossmix_group_velocity.f90
