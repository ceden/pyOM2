


!=======================================================================
! Rayleigh Bernard convection
!=======================================================================

module config_module
 real*8 :: fac=1.0,mix=5e-6 !fac was 2, mix was e-4
 real*8, allocatable :: heat_flux(:,:)
end module config_module

subroutine set_parameter
 use main_module   
 use config_module   
 use diagnostics_module   
 use tke_module   
 implicit none

  nx=int(160*fac) !was 64*
  nz=int(60*fac)
  ny=1
  dt_tracer=0.1!0.01/fac
  dt_mom   =dt_tracer !0.01/fac
     
  enable_conserve_energy = .false.
  coord_degree           = .false.
  enable_cyclic_x        = .true.
  enable_hydrostatic     = .false.
     
  eq_of_state_type       =  1
  enable_tempsalt_sources = .true.

  congr_epsilon = 1e-9
  congr_max_iterations = 5000
  congr_epsilon_non_hydro=   1e-9
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
  enable_diag_snapshots  = .true.; snapint  = 1 !change this for ocean vs.tank
end subroutine set_parameter


subroutine set_grid
 use main_module   
 use config_module   
 implicit none
 dxt(:)=0.005/fac 
 dyt(:)=0.005/fac 
 dzt(:)=0.005/fac 
end subroutine set_grid

subroutine set_coriolis
end subroutine set_coriolis



subroutine set_initial_conditions
 use main_module   
 use config_module   
 implicit none
 integer :: i,k
 !temp(0:160,ny,40:60,:) = 5 !was temp(:,:,:,:)=0
 do k=40,60  
   do i=1,160
    if (i>=is_pe.and.i<=ie_pe)  temp(i,:,k,:) = 5
   enddo
 enddo

 allocate( heat_flux(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx) ); heat_flux = 0
 do i=nx/2-3,nx/2+3
    if (i>=is_pe.and.i<=ie_pe)  heat_flux(i,:) = 5/4185.5/dzt(1)
 enddo

end subroutine set_initial_conditions




subroutine set_forcing
 use main_module
 use config_module   
 implicit none
 !temp_source((nx/2 -3):(nx/2 +3),:,4) =  5/4185.5 /dzt(1) !was temp_source(:,:,1)
 temp_source(:,:,4) = heat_flux
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



