
include ../../site_specific.mk_${HOSTTYPE}

all: main_module.o vertical_friction.o numerics.o diffusion.o bottom_friction.o harmonic_friction.o \
     biharmonic_friction.o momentum.o momentum_advection.o vertical_velocity.o thermodynamics.o restart.o main.o

clean :
	 rm -f *.o *.mod 

main_module.o : main_module.f90 
	$(F90) $(F90FLAGS) -c main_module.f90
numerics.o: numerics.f90 main_module.o 
	$(F90) $(F90FLAGS) -c numerics.f90
vertical_friction.o: vertical_friction.f90 main_module.o 
	$(F90) $(F90FLAGS) -c vertical_friction.f90
bottom_friction.o: bottom_friction.f90 main_module.o 
	$(F90) $(F90FLAGS) -c bottom_friction.f90
harmonic_friction.o: harmonic_friction.f90 main_module.o 
	$(F90) $(F90FLAGS) -c harmonic_friction.f90
biharmonic_friction.o: biharmonic_friction.f90 main_module.o  ../isoneutral/isoneutral_module.o
	$(F90) $(F90FLAGS) -I../isoneutral -c biharmonic_friction.f90
diffusion.o: diffusion.f90 main_module.o 
	$(F90) $(F90FLAGS) -c diffusion.f90

main.o: main.f90 main_module.o ../tke/tke_module.o ../eke/eke_module.o \
        ../idemix/idemix_module.o ../diagnostics/diagnostics_module.o \
        ../etc/timing_module.o
	$(F90) $(F90FLAGS) -I../tke -I../eke -I../idemix -I../diagnostics -I../etc -c main.f90

time_step.o: time_step.f90 main_module.o 
	$(F90) $(F90FLAGS)  $(CDFFLAGS) -c time_step.f90

momentum.o: momentum.f90 main_module.o ../isoneutral/isoneutral_module.o ../etc/timing_module.o \
	     ../idemix/idemix_module.o ../rossmix/rossmix_module.o
	$(F90) $(F90FLAGS) -I../isoneutral -I../idemix -I../rossmix -I../etc -c momentum.f90

momentum_advection.o: momentum_advection.f90 main_module.o 
	$(F90) $(F90FLAGS) -c momentum_advection.f90

vertical_velocity.o: vertical_velocity.f90 main_module.o 
	$(F90) $(F90FLAGS) -c vertical_velocity.f90

thermodynamics.o: thermodynamics.f90 main_module.o \
                  ../isoneutral/isoneutral_module.o ../tke/tke_module.o \
                  ../density/density.o ../etc/timing_module.o ../rossmix/rossmix_module.o
	$(F90) $(F90FLAGS) -I../isoneutral -I../etc -I../density -I../tke -I../rossmix -c thermodynamics.f90

restart.o: restart.f90  main_module.o ../tke/tke_module.o ../eke/eke_module.o \
           ../idemix/idemix_module.o
	$(F90) $(F90FLAGS) -I../tke -I../eke -I../idemix -c restart.f90

