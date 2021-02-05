
include ../../site_specific.mk_${HOSTTYPE}

all: tracer_module.o tracer.o  

tracer_module.o: tracer_module.f90 ../main/main_module.o 
	 	$(F90) $(F90FLAGS) -I../main -c tracer_module.f90
tracer.o: tracer.f90 ../main/main_module.o ../etc/timing_module.o \
               ../isoneutral/isoneutral_module.o ../rossmix/rossmix_module.o tracer_module.o
	 	$(F90) $(F90FLAGS) -I../main -I../isoneutral -I../rossmix -I../etc -c tracer.f90
clean :
	 rm -f *.o *.mod 

