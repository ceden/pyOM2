
!=======================================================================
! 1 D test case for idemix and leewaves
!======================================================================= 

module config_module
 ! use this module only locally in this file
 implicit none
 real*8, parameter :: N_0 = 30e-4,z_b = 800.0  !N_0 = 0.008
 real*8, parameter :: zu0=200.0,zu1=750
 real*8, parameter :: h0=2000.
 real*8, allocatable :: u_target(:,:,:),v_target(:,:,:)
end module config_module


subroutine set_parameter
 ! ----------------------------------
 !       set here main parameter
 ! ----------------------------------
 use main_module   
 use config_module
 use diagnostics_module   
 use idemix_module   

 implicit none

 nx   = 1; nz   = 200; ny   = 1
 coord_degree           = .false.
 enable_cyclic_x        = .true.
 enable_cyclic_y        = .true.
 enable_hydrostatic     = .true.
 eq_of_state_type       = 1
     
 enable_momentum_equation = .true.
 enable_momentum_sources = .true.
 enable_thermodynamic_equation = .false.
 
 runlen =  86400*120
 dt_mom = 360; dt_tracer = dt_mom
 enable_diag_ts_monitor = .true.; ts_monint =  6*3600! dt_mom
 enable_diag_snapshots  = .true.; snapint  =   6*3600! dt_mom
 

 enable_idemix = .true.
 enable_leewaves = .true.
 !enable_idemix_hor_diffusion = .true.
 
 tau_v = 3*86400.
 mu0 = 4./3.
 mu0_min = 0.
 tc_max = 1e-3
 Noverf_min0 = 0.
 
 !idemix3_leewaves_intervall = 2*dt_tracer
 
 enable_leewaves_slowflux = .false.
  
end subroutine set_parameter



subroutine set_grid
 use main_module   
 use config_module   
 implicit none
 dxt    = 5e3
 dyt    = 5e3
 dzt    = h0/nz
end subroutine set_grid

subroutine set_coriolis
 use main_module   
 use config_module
 implicit none
 coriolis_t = 1e-4
end subroutine set_coriolis


subroutine set_initial_conditions
 use main_module
 use tke_module
 use idemix_module   
 use config_module
 implicit none
 integer :: i,j,k
 real*8 :: alpha,get_drhodT,N2


 salt = 35.0
 temp(:,:,nz,:)=25d0
 do k=nz-1,1,-1
   do j=js_pe-onx,je_pe+onx
      do i=is_pe-onx,ie_pe+onx  ! b_z=N^2 , b=-rho/rho_0 g,  T=rho/rho_T, T=-b rho_0/( rho_T g)   
       alpha = get_drhodT(salt(i,j,k,1),temp(i,j,k+1,1),abs(zt(k)))
       N2  = ( 2e-4 + N_0*exp(zw(k)/z_b)  )**2
       !N2  = (  N_0 )**2
       temp(i,j,k,:)  = ( temp(i,j,k+1,:) +dzw(k)*rho_0/(grav*alpha)*N2 )*maskT(i,j,k)
     enddo
   enddo
 enddo


 allocate( u_target(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) )
 allocate( v_target(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) )
 do k=1,nz
   u_target(:,:,k) = 0.1 + 0.1*exp( -(zt(k)+zu1)**2/zu0**2)
   v_target(:,:,k) = 0.0
   u(:,:,k,tau)    =  u_target(:,:,k)
   u(:,:,k,taum1)  =  u_target(:,:,k)
 enddo
 
 h_rms = 10
 k_s = 1./10e3
 ph_s = 0
 k_n = k_s

end subroutine set_initial_conditions


subroutine set_forcing
 ! u_t =   f v - f v*  + Fu_z,  v_t = - f u + f u*  + Fv_z
 use main_module   
 use config_module
 implicit none
 integer :: k
 do k=1,nz
    u_source(:,:,k)     = -coriolis_t(:,:)*v_target(:,:,k)
    v_source(:,:,k)     =  coriolis_t(:,:)*u_target(:,:,k)
 enddo
 !u_source = (.02/86400.)*(u_target-u(:,:,:,tau))
end subroutine set_forcing

subroutine set_topography
 use main_module   
 implicit none
 kbot=1
end subroutine set_topography


subroutine set_diagnostics
end subroutine set_diagnostics

subroutine set_particles
end subroutine set_particles





