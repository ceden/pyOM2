


!=======================================================================
! Rayleigh Bernard convection
!=======================================================================

module config_module
 real*8 :: fac=2.0,mix=5e-3
end module config_module

subroutine set_parameter
 use main_module   
 use config_module   
 use diagnostics_module   
 implicit none

  nx=int(64*fac)
  nz=int(20*fac)
  ny=1
  dt_tracer=0.25/fac
  dt_mom   =0.25/fac
     
  enable_conserve_energy = .false.
  coord_degree           = .false.
  enable_cyclic_x        = .true.
  enable_hydrostatic     = .false.
     
  eq_of_state_type       =  1
  enable_tempsalt_sources = .true.

  congr_epsilon = 1e-6
  congr_max_iterations = 5000
  congr_epsilon_non_hydro=   1e-6
  congr_max_itts_non_hydro = 5000    

  enable_explicit_vert_friction = .true.
  kappam_0 = mix/fac**2
  enable_hor_friction = .true.
  a_h = mix/fac**2

  enable_superbee_advection = .true.
  !enable_hor_diffusion = .true
  !kappah_0 = mix/fac**2
  !k_h = mix/fac**2
     
  runlen =  86400.0
  enable_diag_ts_monitor = .true.; ts_monint =dt_tracer
  enable_diag_snapshots  = .true.; snapint  = 20
end subroutine set_parameter


subroutine set_grid
 use main_module   
 use config_module   
 implicit none
 dxt(:)=0.5/fac 
 dyt(:)=0.5/fac 
 dzt(:)=0.5/fac 
end subroutine set_grid

subroutine set_coriolis
end subroutine set_coriolis



subroutine set_initial_conditions
  use main_module   
  implicit none
  integer :: i,j,k
  do k=1,nz
   do j=js_pe-onx,je_pe+onx
    do i=is_pe-onx,ie_pe+onx
       temp(i,j,k,:)= 0.05*sin(xt(i)/(20*dxt(is_pe))*pi)*maskT(i,j,k)
    enddo
   enddo
  enddo
end subroutine set_initial_conditions




subroutine set_forcing
 use main_module
 implicit none
  temp_source(:,:,nz) = -175/4185.5 /dzt(nz)
  temp_source(:,:,1 ) =  175/4185.5 /dzt(1)
end subroutine set_forcing


subroutine set_topography
 use main_module   
 implicit none
 kbot =1
end subroutine set_topography


subroutine set_diagnostics
end subroutine set_diagnostics


subroutine set_particles
end subroutine set_particles



