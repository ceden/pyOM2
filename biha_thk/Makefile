

include ../../site_specific.mk_${HOSTTYPE}

all: biharmonic_thickness_module.o biharmonic_thickness_mixing.o \
     biharmonic_thickness_friction.o biharmonic_thickness_diag_snap.o \
     biharmonic_thickness_backscatter.o biharmonic_thickness_main.o

clean :
	 rm -f *.o *.mod 


biharmonic_thickness_module.o  : biharmonic_thickness_module.f90 ../main/main_module.o 
	 	$(F90) $(F90FLAGS) -I../main -c biharmonic_thickness_module.f90
        
biharmonic_thickness_main.o: biharmonic_thickness_main.f90 ../main/main_module.o biharmonic_thickness_module.o 
	 	$(F90) $(F90FLAGS)  -I../main  -c biharmonic_thickness_main.f90

biharmonic_thickness_mixing.o: biharmonic_thickness_mixing.f90 ../main/main_module.o \
                                       ../diagnostics/diagnostics_module.o biharmonic_thickness_module.o 
	 	$(F90) $(F90FLAGS)  -I../main -I../diagnostics -c biharmonic_thickness_mixing.f90

biharmonic_thickness_friction.o: biharmonic_thickness_friction.f90 ../main/main_module.o biharmonic_thickness_module.o
	 	$(F90) $(F90FLAGS)  -I../main  -c biharmonic_thickness_friction.f90
        
biharmonic_thickness_backscatter.o: biharmonic_thickness_backscatter.f90 \
                         ../diagnostics/diagnostics_module.o ../main/main_module.o biharmonic_thickness_module.o
	 	$(F90) $(F90FLAGS)  -I../main -I../diagnostics -c biharmonic_thickness_backscatter.f90
        
biharmonic_thickness_diag_snap.o: biharmonic_thickness_diag_snap.f90 ../main/main_module.o   biharmonic_thickness_module.o
	 	$(F90) $(F90FLAGS) $(CDFFLAGS) -I../main  -c biharmonic_thickness_diag_snap.f90


