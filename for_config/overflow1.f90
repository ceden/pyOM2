


!=======================================================================
!
! slope convection
!
!=======================================================================

module config_module
 real*8 :: fac=1.0,mix=5e-3,L0,H0
end module config_module

subroutine set_parameter
 use main_module   
 use config_module   
 use diagnostics_module   
 use tke_module   
 implicit none

  nx=int(128*fac)
  nz=int(40*fac)
  ny=1
  dt_tracer=0.25/fac
  dt_mom   =0.25/fac
     
  enable_conserve_energy = .false.
  coord_degree           = .false.
  enable_hydrostatic     = .false.
     
  eq_of_state_type       =  1

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
 L0 = dxt(is_pe)*nx
 H0 = dzt(1)*nz
end subroutine set_grid

subroutine set_coriolis
end subroutine set_coriolis



subroutine set_initial_conditions
 use main_module   
 use config_module   
 implicit none
  integer :: i
  temp = 0.5
  do i=is_pe-onx,ie_pe+onx
    if (xt(i) <= L0/3) temp(i,:,:,:)=0
  enddo
end subroutine set_initial_conditions




subroutine set_forcing
end subroutine set_forcing


subroutine set_topography
 use main_module   
 use config_module   
 implicit none
 integer :: i,k
 real*8 :: fxa
 kbot =1
 do i=is_pe,ie_pe
    fxa=H0-H0/2*exp( -(xt(i)-L0/3)**2/(L0/4.)**2)
    do k=1,nz
       if (zt(k).lt.-fxa) kbot(i,:)=k
    enddo
 enddo
end subroutine set_topography


subroutine set_diagnostics
end subroutine set_diagnostics


subroutine set_particles
end subroutine set_particles



