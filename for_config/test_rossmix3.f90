
!=======================================================================
! 2 D test case for Rossmix module
!======================================================================= 

module config_module
 ! use this module only locally in this file
 implicit none
 real*8,parameter :: N_0 = 0.012, z_b = 800.0  ! N=N_0 exp(z/z_b)
 real*8, parameter :: U00 = 0.5, y_star=400.e3,dy_star=100e3, h_star=300, dh_star = 500
 real*8, allocatable :: u_ini(:,:,:),v_ini(:,:,:)

 real*8,parameter :: Ek = 0.002
 real*8,parameter :: Ek_res = 0.005
 !real*8,parameter :: Ek = 0.01
 !real*8,parameter :: Ek_res = 0.01
 real*8,parameter :: f0 = 1e-4        ! Coriolis freq.
 real*8,parameter :: dx = 20e3        
 real*8,parameter :: Dz = 2000/20       
 
end module config_module


subroutine set_parameter
 ! ----------------------------------
 !       set here main parameter
 ! ----------------------------------
 use main_module   
 use config_module
 use diagnostics_module   
 use rossmix_module   
 use isoneutral_module   
 use eke_module   
 implicit none
 
 !nx   = 1; nz   = 200; ny   = 1
 nx   = 10; nz   = 20; ny   = 2*20
 dt_mom     = 1200
 dt_tracer  = dt_mom

 coord_degree           = .false.
 enable_cyclic_x        = .true.
 !enable_cyclic_y        = .true.
 enable_hydrostatic     = .true.
 eq_of_state_type       = 1 
 !enable_free_surface = .true.

 !enable_momentum_equation = .false.
 !enable_thermodynamic_equation = .false.
 !enable_conserve_energy = .false. 
 enable_store_lateral_friction_heat = .true.

 congr_epsilon = 1e-12
 congr_max_iterations = 5000
 !enable_streamfunction = .true.
     
 runlen =  365*86400.*10
 enable_diag_snapshots  = .true.; snapint  =  1.0*86400
 enable_diag_ts_monitor = .true.; ts_monint = snapint

 !enable_momentum_sources_with_AB = .true.
 enable_superbee_advection = .true.
 enable_hor_friction = .true.
 A_h      = Ek_res*f0*dx**2 
 !enable_biharmonic_friction    =  .true.
 !a_hbi    = Ek_res*f0*dx**4 
 enable_explicit_vert_friction = .true.
 kappam_0 = Ek_res*f0*dz**2

 !enable_rossmix = .true.
 rossmix_barotropic_fac = 5
 rossmix_calc_cg_int       = 86400.0
 !rossmix_bottom_friction   = 1./(3.0*86400.)
 rossmix_rayleigh_friction   = Ek*f0*10

 ! set of parameter for talk in Copenhagen 2016
 !rossmix_inverse_energy    = 0.2e-12
 !rossmix_barotropisation   = 0.1e-12
 !rossmix_symmetrisation    = 1./(10*86400.)

 ! test new parameter
 rossmix_inverse_energy    = 0.2e-6
 rossmix_barotropisation   = 0.1e-6
 rossmix_symmetrisation    = 1./(20*86400.)

 !enable_rossmix_mean_flow_interaction = .false.
 !enable_rossmix_bolus_form      = .false.
 !enable_rossmix_vertical_stress = .true.
 !enable_rossmix_lateral_stress  = .false.
 enable_rossmix_diffusion = .true.
 rossmix_lat_diffusion = 10.
 rossmix_phi_diffusion =  (2*pi/(np-2))**2 /DX**2*rossmix_lat_diffusion 
 np=20+2
 nm = 2
 ! E_t =  A E_xx ->    dx^2/A = T,   A = dx^2/T


 enable_neutral_diffusion = .true.; 
 !enable_TEM_friction = .true.
 K_iso_0 = 0.0
 K_iso_steep = 1000.0
 iso_dslope=4./1000.0
 iso_slopec=4./1000.0
 K_gm_0 = 500.0

 !enable_eke = .true.
 !eke_c_eps  = 0.5
 !eke_k_max  = 1e4
 !eke_c_k    = 0.4
 !eke_c_eps  = 0.5
 !eke_cross  = 2.
 !eke_crhin  = 1.0
 !eke_lmin   = 100.0
 !enable_eke_superbee_advection  = .true.
 !enable_eke_isopycnal_diffusion = .true.
 !enable_eke_diss_surfbot = .true.
 !eke_diss_surfbot_frac = 0.2

end subroutine set_parameter



subroutine set_grid
 use main_module   
 use config_module   
 implicit none
 dxt    = Dx
 dyt    = Dx
 dzt    = Dz !H0/nz !10.0
 y_origin= Dx
end subroutine set_grid


subroutine set_coriolis
  use main_module   
  use config_module   
  implicit none
  integer :: j
  do j=js_pe-onx,je_pe+onx
    coriolis_t(:,j) = f0 + 2d-11*yt(j)
  enddo
end subroutine set_coriolis



subroutine set_initial_conditions
 use main_module
 use rossmix_module   
 use eke_module   
 use config_module
 implicit none
 integer :: i,j,k
 real*8 :: alpha,get_drhodT,N2
 real*8 :: aloc(nx,ny),bloc(ny,nz),cloc(ny,nz)
 real*8 :: uz(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz)
! integer,parameter :: seed = 86456
 integer,parameter :: seed = 86457
 real*4 :: rand

 alpha = get_drhodT(35d0,5d0,0d0) 

 allocate( u_ini(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); u_ini=0.
 allocate( v_ini(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); v_ini=0.
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

 do j=js_pe,je_pe
  do k=1,nz
   u_ini(:,j,k)     = U00*exp(-(yt(j)-y_star)**2/dy_star**2)*0.5*(1+tanh( (zt(k)+h_star)/dh_star ))
   !uz   (:,j,k)     = U00/dh_star*exp(-(yt(j)-y_star)**2/dy_star**2)*(1-tanh( (zt(k)+h_star)/dh_star )**2)
  enddo
 enddo
 do k=2,nz-1
   uz(:,:,k) = (u_ini(:,:,k+1)-u_ini(:,:,k-1))/(dzt(k+1)+dzt(k))
 enddo
 uz(:,:,1)=uz(:,:,2)
 uz(:,:,nz)=uz(:,:,nz-1)


 u(:,:,:,tau)   = u_ini; !v(:,:,:,tau)   = v_ini
 u(:,:,:,taum1) = u_ini; !v(:,:,:,taum1) = v_ini
 u(:,:,:,taup1) = u_ini; !v(:,:,:,taup1) = v_ini

 ! rho0 fu = -p_y, p_z = -g rho, rho0 f u_z = -g rho_y,  rho_y = - rho0 f u_z/g = alpha T_y
 do k=1,nz
   aloc(1,js_pe:je_pe) = dyt(js_pe:je_pe)*uz(is_pe,js_pe:je_pe,k) *coriolis_t(is_pe,js_pe:je_pe)/grav/alpha*rho_0
   call pe0_recv_2D(nx,ny,aloc)
   cloc(1:ny,k) = aloc(1,1:ny)
 enddo
 bloc(1,:)=0.
 do j=1,ny-1
    bloc(j+1,:)=bloc(j,:) + cloc(j,:)
    call pe0_bcast(bloc(j+1,:),nz)
 enddo

 call srand(seed)

 do j=js_pe,je_pe
   do k=1,nz
    do i=is_pe,ie_pe
      temp(i,j,k,:)  = temp(i,j,k,:) + bloc(j,k) + 1e-3*rand()
    enddo
   enddo
 enddo


 if (enable_rossmix) then
  do k=2,np-1
    do j=js_pe,je_pe
     do i=is_pe,ie_pe
      E_s(i,j,k,:,1:nm)= 1e-3*rand()
      E_l(i,j,k,:,1:nm)= 1e-3*rand()
      !E_l(i,j,k,:,1)= 1e3*exp( -(yt(j)-500e3)**2/50e3**2  &
      !                         -( phit(k)-3*pi/4.)**2/ (pi/5.)**2 ) *maskTp(i,j,k) 
     enddo
    enddo
  enddo
 endif

 if (enable_eke) then
  do k=1,nz
    do j=js_pe,je_pe
     do i=is_pe,ie_pe
      eke(i,j,k,:)= 1e-8*rand()
     enddo
    enddo
  enddo
 endif

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





