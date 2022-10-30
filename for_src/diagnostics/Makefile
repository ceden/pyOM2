
include ../../site_specific.mk_${HOSTTYPE}

all: diagnostics_module.o diag_main.o diag_averages.o diag_over.o diag_snap.o \
     diag_energy.o diag_particles.o diag_particles_netcdf.o diag_variance.o \
     diag_snap_tke.o diag_snap_eke.o diag_snap_isoneutral.o diag_snap_energy.o \
     diag_snap_idemix.o diag_snap_tracer.o diag_snap_rossmix.o diag_snap_fluxes.o \
     diag_snap_tendency.o diag_snap_parallel.o diag_snap_chunks.o

clean :
	 rm -f *.o *.mod 

diagnostics_module.o: diagnostics_module.f90
	$(F90) $(F90FLAGS) -c diagnostics_module.f90

diag_main.o: diag_main.f90 ../main/main_module.o ../tke/tke_module.o \
             ../isoneutral/isoneutral_module.o ../eke/eke_module.o ../rossmix/rossmix_module.o\
             ../idemix/idemix_module.o ../tracer/tracer_module.o diagnostics_module.o
	$(F90) $(F90FLAGS) $(CDFFLAGS) -I../main -I../tke -I../isoneutral -I../eke -I../idemix -I../rossmix -I../tracer -c diag_main.f90

diag_averages.o: diag_averages.f90 ../main/main_module.o diagnostics_module.o
	$(F90) $(F90FLAGS) $(CDFFLAGS) -I../main -c diag_averages.f90

diag_variance.o: diag_variance.f90 ../main/main_module.o diagnostics_module.o
	$(F90) $(F90FLAGS) $(CDFFLAGS) -I../main -c diag_variance.f90

diag_energy.o: diag_energy.f90 ../main/main_module.o ../tke/tke_module.o ../isoneutral/isoneutral_module.o \
               ../eke/eke_module.o ../idemix/idemix_module.o diagnostics_module.o
	$(F90) $(F90FLAGS) $(CDFFLAGS) -I../main -I../tke -I../eke -I../idemix -I../isoneutral -c diag_energy.f90

diag_snap.o: diag_snap.f90 ../main/main_module.o diagnostics_module.o
	$(F90) $(F90FLAGS) $(CDFFLAGS) -I../main  -c diag_snap.f90

diag_snap_idemix.o: diag_snap_idemix.f90 ../main/main_module.o ../idemix/idemix_module.o diagnostics_module.o
	$(F90) $(F90FLAGS) $(CDFFLAGS) -I../main  -I../idemix  -c diag_snap_idemix.f90

diag_snap_rossmix.o: diag_snap_rossmix.f90 ../main/main_module.o ../rossmix/rossmix_module.o diagnostics_module.o
	$(F90) $(F90FLAGS) $(CDFFLAGS) -I../main  -I../rossmix  -c diag_snap_rossmix.f90

diag_snap_energy.o: diag_snap_energy.f90 ../main/main_module.o ../isoneutral/isoneutral_module.o ../idemix/idemix_module.o  diagnostics_module.o
	$(F90) $(F90FLAGS) $(CDFFLAGS) -I../main  -I../isoneutral  -I../idemix  -c diag_snap_energy.f90

diag_snap_isoneutral.o: diag_snap_isoneutral.f90 ../main/main_module.o ../isoneutral/isoneutral_module.o  diagnostics_module.o
	$(F90) $(F90FLAGS) $(CDFFLAGS) -I../main -I../isoneutral -c diag_snap_isoneutral.f90

diag_snap_eke.o: diag_snap_eke.f90 ../main/main_module.o ../eke/eke_module.o diagnostics_module.o
	$(F90) $(F90FLAGS) $(CDFFLAGS) -I../main  -I../eke -c diag_snap_eke.f90

diag_snap_tracer.o: diag_snap_tracer.f90 ../main/main_module.o ../tracer/tracer_module.o diagnostics_module.o
	$(F90) $(F90FLAGS) $(CDFFLAGS) -I../main  -I../tracer  -c diag_snap_tracer.f90

diag_snap_tke.o: diag_snap_tke.f90 ../main/main_module.o ../tke/tke_module.o diagnostics_module.o
	$(F90) $(F90FLAGS) $(CDFFLAGS) -I../main -I../tke -c diag_snap_tke.f90

diag_over.o: diag_over.f90 ../main/main_module.o ../isoneutral/isoneutral_module.o 
	$(F90) $(F90FLAGS) $(CDFFLAGS) -I../main -I../isoneutral -c diag_over.f90

diag_particles.o: diag_particles.f90 ../main/main_module.o 
	$(F90) $(F90FLAGS) -I../main -c diag_particles.f90
diag_particles_netcdf.o: diag_particles_netcdf.f90 ../main/main_module.o  diag_particles.o
	$(F90) $(F90FLAGS) $(CDFFLAGS) -I../main -c diag_particles_netcdf.f90

diag_snap_fluxes.o: diag_snap_fluxes.f90 ../main/main_module.o  ../isoneutral/isoneutral_module.o 
	$(F90) $(F90FLAGS) $(CDFFLAGS) -I../main  -I../isoneutral  -c diag_snap_fluxes.f90

diag_snap_tendency.o: diag_snap_tendency.f90 ../main/main_module.o  ../isoneutral/isoneutral_module.o 
	$(F90) $(F90FLAGS) $(CDFFLAGS) -I../main  -I../isoneutral  -c diag_snap_tendency.f90
    
diag_snap_parallel.o: diag_snap_parallel.F90  ../main/main_module.o
	$(F90) $(F90FLAGS) $(CDFFLAGS)  $(MPIFLAGS)  -I../main -c diag_snap_parallel.F90
    
diag_snap_chunks.o: diag_snap_chunks.f90 ../main/main_module.o
	$(F90) $(F90FLAGS) $(CDFFLAGS)   -I../main -c diag_snap_chunks.f90
    
