
  
!=======================================================================
! Internal wave maker
!=======================================================================

module config_module
 use main_module   
 real*8 :: fac=2.0
 real*8 :: N_0 = 2*pi/10.
 real*8 :: OM0 = 1./(1.5*10),x0
 real*8, allocatable :: t0(:,:,:),dt0(:,:,:,:),u0(:,:,:)
end module config_module

subroutine set_parameter
 use main_module   
 use config_module   
 use diagnostics_module   
 use tke_module   
 implicit none
  ny=1
  nx=int(64*fac)
  nz=int(64*fac)

  dt_mom   = 20*0.025/fac
  dt_tracer= 20*0.025/fac
     
  enable_conserve_energy = .false.
  coord_degree           = .false.
  enable_cyclic_x        = .true.
  enable_hydrostatic     = .false.
  eq_of_state_type       = 1

  congr_epsilon = 1e-12
  congr_max_iterations = 5000
  congr_epsilon_non_hydro=   1e-9
  congr_max_itts_non_hydro = 5000    

  enable_explicit_vert_friction = .true.
  kappam_0 = 5e-3/fac**2   
  enable_hor_friction =  .true.
  a_h      = 5e-3/fac**2
  enable_superbee_advection = .true.

  enable_tempsalt_sources = .true.
  enable_momentum_sources = .true.

  runlen =  86400.0
  enable_diag_ts_monitor = .true.; ts_monint = 0.5!dt_mom
  enable_diag_snapshots  = .true.; snapint  = 5.0!dt_mom

end subroutine set_parameter


subroutine set_grid
 use main_module   
 use config_module   
 implicit none
 dxt(:)=0.25/fac 
 dyt(:)=0.25/fac 
 dzt(:)=0.25/fac 
end subroutine set_grid



subroutine set_initial_conditions
  use main_module   
  use config_module   
  implicit none
  integer :: i,k
  real*8 :: alpha,get_drhodT

  alpha = get_drhodT(35d0,5d0,0d0) 
  allocate( t0(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); t0 = 0
  allocate( u0(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); u0 = 0
  allocate( dt0(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz,3) ); dt0 = 0
  do k=1,nz
    t0(:,:,k)=-N_0**2*zt(k)/9.81/alpha*rho_0*maskT(:,:,k)
  enddo
  x0 = 0.0; if (nx/2 >= is_pe .and. nx/2<= ie_pe) x0 = xu(nx/2)
  call global_max(x0)
  do k=1,nz
   do i=is_pe,ie_pe
    u0(i,:,k)= maskU(i,:,k)*1./(100*60.*dt_tracer)*exp( -(xu(i)-x0)**2/(dxt(is_pe)*1)**2 ) &
                                                  *exp( -(zt(k)-zt(nz/2))**2/(dzt(1)*1)**2 )
   enddo
  enddo
end subroutine set_initial_conditions




subroutine set_forcing
 use main_module
  use config_module   
 implicit none

 ! implement effect of background state
 
 ! update density, etc of last time step
 temp(:,:,:,tau) = temp(:,:,:,tau) + t0
 call calc_eq_of_state(tau)
 temp(:,:,:,tau) = temp(:,:,:,tau) - t0
      
 ! advection of background temperature
 call advect_tracer(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,t0,dt0(:,:,:,tau)) 
 temp_source = (1.5+ab_eps)*dt0(:,:,:,tau) - ( 0.5+ab_eps)*dt0(:,:,:,taum1)

 ! wavemaker
 u_source= u0*sin(2*pi*OM0*itt*dt_tracer)
 
end subroutine set_forcing

subroutine set_coriolis
end subroutine set_coriolis

subroutine set_topography
 use main_module   
 implicit none
 kbot =1
end subroutine set_topography


subroutine set_diagnostics
end subroutine set_diagnostics


subroutine set_particles
end subroutine set_particles



