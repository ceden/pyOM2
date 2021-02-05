
include ../../site_specific.mk_${HOSTTYPE}

all: adv_main.o adv_2nd.o adv_superbee.o adv_upwind.o adv_upwind3.o adv_dst3.o

clean :
	 rm -f *.o *.mod 

adv_main.o: adv_main.f90 ../main/main_module.o 
	$(F90) $(F90FLAGS) -I../main -c adv_main.f90

adv_2nd.o: adv_2nd.f90 ../main/main_module.o 
	$(F90) $(F90FLAGS) -I../main -c adv_2nd.f90

adv_superbee.o: adv_superbee.f90 ../main/main_module.o 
	$(F90) $(F90FLAGS) -I../main -c adv_superbee.f90

adv_upwind.o: adv_upwind.f90 ../main/main_module.o 
	$(F90) $(F90FLAGS) -I../main -c adv_upwind.f90

adv_upwind3.o: adv_upwind3.f90 ../main/main_module.o 
	$(F90) $(F90FLAGS) -I../main -c adv_upwind3.f90

adv_dst3.o: adv_dst3.f90 ../main/main_module.o 
	$(F90) $(F90FLAGS) -I../main -c adv_dst3.f90

