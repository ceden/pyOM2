
include ../../site_specific.mk_${HOSTTYPE}

all: obc_module.o obc.o

obc_module.o: obc_module.f90 ../main/main_module.o 
	$(F90) $(F90FLAGS) -I../main -c obc_module.f90
obc.o: obc.f90 ../main/main_module.o obc_module.o 
	$(F90) $(F90FLAGS) -I../main -c obc.f90
clean :
	 rm -f *.o *.mod 

