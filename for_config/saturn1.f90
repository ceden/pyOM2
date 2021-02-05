
!=======================================================================
!  
!======================================================================= 

module config_module
 implicit none
 real*8, allocatable :: U_0(:,:,:),T_0(:,:,:)
 real*8,parameter :: Ek  = 0.01
 real*8,parameter :: deg_wide = 20.
 real*8,parameter :: delta_x = 0.25
 real*8,parameter :: H0 = 2000.0
 real*8,parameter :: N0 = 40*1e-4
 real*8,parameter :: U0 = 0.5
 real*8,parameter :: M0 = sqrt(1e-4*U0/H0)
 real*8,parameter :: CFL = 0.01
end module config_module


subroutine set_parameter
 ! ----------------------------------
 !       set here main parameter
 ! ----------------------------------
 use main_module   
 use diagnostics_module  
 use config_module 
 implicit none
 
  nx   = int(deg_wide/delta_x*2); nz   = 10.; ny  = nx
  
  coord_degree  = .true. 
  congr_epsilon = 1e-12
  congr_max_iterations = 5000
  enable_free_surface = .true.
     
  enable_hydrostatic      = .true.
  enable_cyclic_x         = .true.
  enable_cyclic_y         = .true.
     
  enable_superbee_advection     = .true.
  enable_explicit_vert_friction = .true.
  enable_biharmonic_friction    =  .true.

  enable_conserve_energy =  .false.
  
  runlen =  365*86400.*100

  enable_diag_ts_monitor = .true.; ts_monint = 86400.
  enable_diag_snapshots  = .true.; snapint  =  86400.
  
  eq_of_state_type = 1 
  
  enable_tempsalt_sources = .true.
 
end subroutine set_parameter


subroutine set_grid
  use main_module   
   use config_module
  implicit none
  dzt = H0/nz
  dxt = delta_x
  dyt = delta_x
  x_origin= -deg_wide
  y_origin= -deg_wide
  
  a_hbi    = Ek*2*omega*dxt(is_pe)**4 
  kappam_0 = Ek*2*omega*dzt(1)**2
  
  dt_mom    = CFL/U0*dxt(is_pe)*degtom  
  dt_tracer = dt_mom 
end subroutine set_grid


subroutine set_coriolis
 use main_module   
 implicit none
 integer :: i,j
 do j=js_pe-onx,je_pe+onx
  do i=is_pe-onx,ie_pe+onx
   coriolis_t(i,j) = 2*omega*cos( yt(j)/180.*pi )*cos( xt(i)/180.*pi ) 
  enddo 
 enddo
end subroutine set_coriolis


subroutine set_initial_conditions
 use main_module
 use config_module


 implicit none
 
 integer :: k,j,i
 real*8 :: alpha,get_drhodT,fxa
 real*8 :: Ly
 
 if (my_pe==0) then
 print*,' Rossby radius ', N0*H0/coriolis_t(is_pe,js_pe) /1e3,' km'
 print*,' delta x       ', dxt(is_pe)/1e3,'km'
 print*,' delta t       ', dt_mom,'s'
 endif
 
 alpha = get_drhodT(35d0,5d0,0d0) 

 allocate( U_0(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) )
 allocate( T_0(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) )
 U_0 = 0; T_0 = 0
    
 ! zonal velocity 
 do k=1,nz
    U_0(:,:,k)= 0. 
    u(:,:,k,tau)= U_0(:,:,k)
 enddo     
 u(:,:,:,taum1) = u(:,:,:,tau)

 
 do k=1,nz
     T_0(:,:,k)=-N0**2*zt(k)/grav/alpha*rho_0*maskT(:,:,k)
 enddo   
 
 Ly = nx*dyt(js_pe) 
 fxa = M0**2*Ly/grav/alpha*rho_0
 do j=js_pe,je_pe
   do i=is_pe,ie_pe
    T_0(i,j,:) = T_0(i,j,:) - fxa*exp(-(yt(j)**2+xt(i)**2)/(0.5*deg_wide)**2)
   enddo 
 enddo

 do k=1,nz
  temp(:,:,k,tau)   = T_0(:,:,k)
  temp(:,:,k,taum1) = T_0(:,:,k)
 enddo

end subroutine set_initial_conditions



subroutine set_forcing

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



