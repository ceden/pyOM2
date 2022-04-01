
include ../../site_specific.mk_${HOSTTYPE}


convert_restart.o:  convert_restart.f90
	$(F90) $(F90FLAGS) convert_restart.f90 -o convert_restart.x

clean :
	 rm -f *.o *.mod 

