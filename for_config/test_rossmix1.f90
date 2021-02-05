
! idealized setup to test wave propagation in Rossmix

module config_module
 ! use this module only locally in this file
 implicit none
 real*8,parameter :: N_0 = 0.012, z_b = 800.0  ! N=N_0 exp(z/z_b)

 real*8, parameter :: L0 = 1000e3
 real*8, parameter :: H0 = 2000
 real*8, parameter :: f0 = 0.5e-4      
 real*8, parameter :: BETA0 = 2e-11   

 real*8, parameter :: DX   = L0/32
 real*8, parameter :: DZ   = H0/10
 real*8, parameter :: DT   = 1200.0
 
end module config_module





subroutine set_parameter
 use main_module   
 use diagnostics_module   
 use rossmix_module   
 use config_module
 implicit none

  nx=int(L0/DX); ny=int(L0/DX); nz = int(H0/DZ)
  dt_tracer= DT
  dt_mom = dt_tracer

  runlen=86400.*365.*1000
  enable_diag_ts_monitor = .true.; ts_monint = dt_tracer!*20
  enable_diag_snapshots  = .true.; snapint  =  8640

  enable_cyclic_x  = .false.
  coord_degree     = .false.
  eq_of_state_type = 1
  enable_momentum_equation = .false.
  enable_thermodynamic_equation = .false.

  enable_rossmix = .true.
  rossmix_barotropic_fac =10 
  enable_rossmix_mean_flow_interaction = .false.
  enable_rossmix_bolus_form      = .false.
  enable_rossmix_vertical_stress = .false.
  enable_rossmix_lateral_stress = .false.
  enable_rossmix_diffusion = .false.
  rossmix_calc_cg_int       = 10*86400.0
  np = 40
  nm = 2
end subroutine set_parameter


subroutine set_grid
 use main_module   
 use rossmix_module   
 use config_module
 implicit none
 dxt    = DX
 dyt    = DX
 dzt    = DZ
 y_origin= DX
end subroutine set_grid


subroutine set_coriolis
 use main_module   
 use rossmix_module   
 use config_module
 implicit none
 integer :: j
 do j=js_pe-onx,je_pe+onx
    coriolis_t(:,j) = f0 + BETA0*yt(j)
 enddo
end subroutine set_coriolis

 
subroutine set_topography
  use main_module   
  implicit none
  kbot=1
end subroutine set_topography


subroutine set_initial_conditions
 use main_module
 use rossmix_module
 use config_module
 implicit none
 integer :: i,j,k
 real*8 :: alpha,get_drhodT,N2
 integer,parameter :: seed = 86456


 alpha = get_drhodT(35d0,5d0,0d0) 
 salt = 35.0
 temp(:,:,nz,:)=25d0
 do k=nz-1,1,-1
   do j=js_pe-onx,je_pe+onx
      do i=is_pe-onx,ie_pe+onx  ! b_z=N^2 , b=-rho/rho_0 g,  T=rho/rho_T, T=-b rho_0/( rho_T g)   
       N2 = ( N_0*exp(zt(k)/z_b)  )**2
       temp(i,j,k,:)  = ( temp(i,j,k+1,:) +dzw(k)*rho_0/(grav*alpha)*N2 )*maskT(i,j,k)
     enddo
   enddo
 enddo

 do k=2,np-1
  do j=js_pe,je_pe
    do i=is_pe,ie_pe
      E_s(i,j,k,:,:)= 10*exp( -( xt(i)-500e3)**2/50e3**2 - (yt(j)-500e3)**2/50e3**2  &
                               -( phit(k)-4*pi/4.)**2/(pi/10)**2 ) *maskTp(i,j,k) 
    enddo
  enddo
 enddo
end subroutine set_initial_conditions


subroutine set_forcing
end subroutine set_forcing

subroutine set_diagnostics
end subroutine set_diagnostics

subroutine set_particles
end subroutine set_particles
