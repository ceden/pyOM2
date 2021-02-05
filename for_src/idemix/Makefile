
include ../../site_specific.mk_${HOSTTYPE}

all: idemix_module.o idemix.o idemix_low_mode.o idemix3.o\
     idemix_spectral.o idemix_group_vel.o idemix_friction.o idemix_leewaves.o

clean :
	 rm -f *.o *.mod 

idemix_module.o  : idemix_module.f90 ../main/main_module.o 
	 	$(F90) $(F90FLAGS) -I../main -c idemix_module.f90
idemix.o: idemix.f90 ../main/main_module.o ../eke/eke_module.o ../etc/timing_module.o idemix_module.o
	 	$(F90) $(F90FLAGS) -I../main -I../eke -I../etc -c idemix.f90
idemix_low_mode.o: idemix_low_mode.f90 ../main/main_module.o idemix_module.o
	 	$(F90) $(F90FLAGS) -I../main -c idemix_low_mode.f90
idemix_spectral.o: idemix_spectral.f90 ../main/main_module.o idemix_module.o
	 	$(F90) $(F90FLAGS) -I../main -c idemix_spectral.f90
idemix_group_vel.o: idemix_group_vel.f90 ../main/main_module.o idemix_module.o
	 	$(F90) $(F90FLAGS) -I../main -c idemix_group_vel.f90
idemix3.o : idemix3.f90 ../main/main_module.o  idemix_module.o
		$(F90) $(F90FLAGS) -I../main -c idemix3.f90
idemix_friction.o : idemix_friction.f90 ../main/main_module.o idemix_module.o
		$(F90) $(F90FLAGS) -I../main -c idemix_friction.f90

idemix_leewaves.o : idemix_leewaves.f90 ../main/main_module.o idemix_module.o
		$(F90) $(F90FLAGS) -I../main  -c idemix_leewaves.f90
