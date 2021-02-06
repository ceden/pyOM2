


include ../site_specific.mk_$(HOSTTYPE)

code = ../for_src
main = $(code)/main
external = $(code)/external
iso = $(code)/isoneutral
ide = $(code)/idemix
adv = $(code)/advection
dia = $(code)/diagnostics
rossmix = $(code)/rossmix


for_files = $(main)/main_module.f90  \
            $(ide)/idemix_module.f90 \
            $(code)/density/density.f90 \
            $(code)/tke/tke_module.f90 \
            $(code)/eke/eke_module.f90 \
            $(code)/obc/obc_module.f90 \
            $(code)/tracer/tracer_module.f90 \
            $(iso)/isoneutral_module.f90  \
            $(code)/etc/timing_module.f90  \
            $(rossmix)/rossmix_module.f90 \
            $(main)/momentum.f90 $(main)/biharmonic_friction.f90 $(main)/numerics.f90 \
            $(main)/harmonic_friction.f90 $(main)/bottom_friction.f90  $(main)/vertical_friction.f90 \
            $(main)/thermodynamics.f90 $(main)/momentum_advection.f90 $(main)/vertical_velocity.f90 \
            $(adv)/adv_main.f90 \
            $(adv)/adv_2nd.f90 $(adv)/adv_dst3.f90 \
            $(adv)/adv_superbee.f90 $(adv)/adv_upwind.f90 $(adv)/adv_upwind3.f90 \
            $(main)/diffusion.f90  $(main)/restart.f90 \
            $(code)/non_hydrostatic/non_hydrostatic.f90 \
            $(code)/obc/obc.f90  $(code)/density/get_rho.f90  \
            $(external)/solve_pressure.f90 $(external)/streamfct_init.f90 \
            $(external)/solve_streamfct.f90 $(external)/island_module.f90 \
            $(iso)/isoneutral_diffusion.f90 $(iso)/isoneutral_friction.f90 \
            $(iso)/isoneutral.f90 $(iso)/isoneutral_ml.f90 \
            $(ide)/idemix.f90 $(ide)/idemix_low_mode.f90 $(ide)/idemix_friction.f90 \
            $(ide)/idemix_spectral.f90 $(ide)/idemix_group_vel.f90 $(ide)/idemix3.f90 $(ide)/idemix_leewaves.f90 \
            $(rossmix)/rossmix.f90 $(rossmix)/rossmix_flux.f90 $(rossmix)/rossmix_integrate.f90  \
            $(rossmix)/rossmix_friction.f90 $(rossmix)/rossmix_group_velocity.f90  \
            $(rossmix)/rossmix_reflect.f90  \
			$(code)/tke/tke.f90 $(code)/eke/eke.f90 $(code)/tracer/tracer.f90 \
            $(dia)/diag_particles.f90 $(code)/etc/py_stuff.f90 


all : pyOM_code.so pyOM_code_MPI.so

pyOM_code.so: $(for_files) $(code)/parallel/parallel_none.f90
	f2py -c  $(for_files) $(code)/parallel/parallel_none.f90 -m pyOM_code ${F2PYFLAGS}

pyOM_code_MPI.so: $(for_files) $(code)/parallel/parallel_mpi.f90
	f2py -c  $(for_files) $(code)/parallel/parallel_mpi.f90 -m pyOM_code_MPI ${F2PYFLAGS} ${F2PY_MPIFLAGS}

clean: 
	(cd ${code}; make clean)
	rm -f *.so


