
!=======================================================================
! 2 D test case for Rossmix module
!======================================================================= 

module config_module
 ! use this module only locally in this file
 implicit none
 real*8,parameter :: N_0 = 0.012, z_b = 800.0  ! N=N_0 exp(z/z_b)
 real*8, parameter :: U00 = 0.5, y_star=400.e3,dy_star=100e3, h_star=300, dh_star = 500
 real*8, allocatable :: u_ini(:,:,:),v_ini(:,:,:)

 real*8,parameter :: BETA0 = 2e-11

 real*8,parameter :: hfac = 2.0
 real*8,parameter :: vfac = 4.0
 real*8,parameter :: Ek = 0.001
 real*8,parameter :: f0 = 1e-4        ! Coriolis freq.
 real*8,parameter :: dx = 40e3 /hfac
 real*8,parameter :: Dz = 2000/(10*vfac)

 ! Ek = r/f   r = Ek*f
 real*8,parameter :: u_rest = Ek*f0   ! restoring time scale
 
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
 
 nx   = 1; nz   = int(10*vfac); ny   = int(20*hfac)
 dt_mom     = 1800.0/hfac
 dt_tracer  = dt_mom

 coord_degree           = .false.
 enable_cyclic_x        = .true.
 enable_hydrostatic     = .true.
 eq_of_state_type       = 1 
 enable_momentum_sources = .true.
 !enable_free_surface = .true.

 !enable_momentum_equation = .false.
 !enable_thermodynamic_equation = .false.
 !enable_conserve_energy = .false. 
 enable_store_lateral_friction_heat = .true.

 congr_epsilon = 1e-12
 congr_max_iterations = 5000
 !enable_streamfunction = .true.
    
 runlen =  365*86400.*10
 enable_diag_snapshots  = .true.; snapint  =  86400
 enable_diag_ts_monitor = .true.; ts_monint = snapint

 !enable_momentum_sources_with_AB = .true.
 enable_superbee_advection = .true.
 enable_hor_friction = .true.
 A_h      = Ek*f0*dx**2 
 !enable_biharmonic_friction    =  .true.
 !a_hbi    = Ek*f0*dx**4 
 enable_explicit_vert_friction = .true.
 kappam_0 = Ek*f0*dz**2

 enable_rossmix = .true.

 ! max group velocity is beta*R0**2 = beta*(9.81*H0/f0**2)
 rossmix_barotropic_fac = ceiling(  dt_mom* BETA0*grav*(nz*Dz)/f0**2 /DX*2 )

 rossmix_calc_cg_int       = 86400.0
 rossmix_rayleigh_friction   = Ek*f0
 rossmix_N2min = 5e-6
 rossmix_inverse_energy    = 0.5e-6
 rossmix_barotropisation   = 3.0e-6
 !rossmix_symmetrisation    = 1./(20*86400.)
 rossmix_c_tau    = 0.2
 rossmix_c_lambda = 6.0
 enable_rossmix_use_gamma = .true.


 !enable_rossmix_mean_flow_interaction = .false.
 !enable_rossmix_bolus_form      = .false.
 !enable_rossmix_vertical_stress = .true.
 enable_rossmix_lateral_stress  = .false.
 enable_rossmix_diffusion = .true.
 rossmix_lat_diffusion = 25.
 rossmix_phi_diffusion =  (2*pi/(nphi-2))**2 /DX**2*rossmix_lat_diffusion 
 nphi=20+2
 nmodes = 2
 ! E_t =  A E_xx ->    dx^2/A = T,   A = dx^2/T


 !enable_neutral_diffusion = .true.; 
 !!!enable_TEM_friction = .true.
 !enable_skew_diffusion = .true.
 K_iso_0 = 0.0
 K_iso_steep = 000.0
 iso_dslope=1/1000.0
 iso_slopec=10/1000.0
 K_gm_0 = 2500.0

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
    coriolis_t(:,j) = f0 + BETA0*yt(j)
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
 ! real*8 :: aloc(nx,ny),bloc(ny,nz),cloc(ny,nz)
 !real*8 :: dloc(ny,nz),bloc(ny,nz),cloc(ny,nz)
 real*8 :: aloc(is_pe-onx:ie_pe+onx,1:ny,nz)
 real*8 :: bloc(is_pe-onx:ie_pe+onx,1:ny,nz)
 real*8 :: cloc(nx,ny)
 real*8 :: uz(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz)
 integer,parameter :: seed = 86456
! integer,parameter :: seed = 86457
 real*4 :: rand
 !integer :: seed
 !namelist /nml_seed/ seed

 
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

 do j=js_pe-onx,je_pe+onx
  do k=1,nz
   u_ini(:,j,k)     = U00*exp(-(yt(j)-y_star)**2/dy_star**2)*0.5*(1+tanh( (zt(k)+h_star)/dh_star ))
  enddo
 enddo

 uz(:,:,1) = 0.
 do k=1,nz
  uz(:,:,1) = uz(:,:,1) + u_ini(:,:,k)
 enddo
 do k=1,nz
  u_ini(:,:,k) = u_ini(:,:,k) - uz(:,:,1)/nz
 enddo

 u(:,:,:,tau)   = u_ini; !v(:,:,:,tau)   = v_ini
 u(:,:,:,taum1) = u_ini; !v(:,:,:,taum1) = v_ini
 u(:,:,:,taup1) = u_ini; !v(:,:,:,taup1) = v_ini

 ! rho0 fu = -p_y, p_z = -g rho, rho0 f u_z = -g rho_y,  rho_y = - rho0 f u_z/g = alpha T_y

 ! discrete for equidistant grid
 !
 ! (p(k+1)-p(k))/dz  =  -(rho(k+1)+rho(k))*grav/rho_0/2
 !
 ! (f(j)*(u(i-1,j,k)+u(i,j,k)) + f(j+1)*(u(i-1,j+1,k)+u(i,j+1,k)) )/4  = - (p(i,j+1,k)-p(i,j,k))/dy
 !
 ! with uz(k) = (u(k+1)-u(k))/dz
 !
 ! (f(j)*(uz(i-1,j,k)+uz(i,j,k)) + f(j+1)*(uz(i-1,j+1,k)+uz(i,j+1,k)) )*dy*rho_0/grav/2 
 !  +rho(i,j,k+1) - rho(i,j+1,k+1) + rho(i,j,k)  = rho(i,j+1,k) 


 do k=1,nz-1
   uz(:,:,k) = (u_ini(:,:,k+1)-u_ini(:,:,k))/dzw(k)
 enddo
 do k=1,nz-1
  do j=js_pe,je_pe
   do i=is_pe,ie_pe
    cloc(i,j) =dyt(j)*rho_0/grav/2./alpha*( coriolis_t(i,j  )*(uz(i-1,j  ,k)+uz(i,j,k)) &
                                          + coriolis_t(i,j+1)*(uz(i-1,j+1,k)+uz(i,j+1,k)) )
   enddo
  enddo
  call pe0_recv_2D(nx,ny,cloc)
  call pe0_bcast(cloc,nx*ny)
  aloc(is_pe:ie_pe,1:ny,k) = cloc(is_pe:ie_pe,1:ny)
 enddo

 bloc=0
 do k=1,nz-1
  do j=1,ny-1
   bloc(:,j+1,k+1) = bloc(:,j,k) + bloc(:,j,k+1) + aloc(:,j,k) - bloc(:,j+1,k) 
  enddo
 enddo

 do k=1,nz
  do j=js_pe,je_pe
    temp(:,j,k,tau)   = temp(:,j,k,tau)   + bloc(:,j,k)
    temp(:,j,k,taum1) = temp(:,j,k,taum1) + bloc(:,j,k)
  enddo
 enddo

 !do k=2,nz-1
 !  uz(:,:,k) = (u_ini(:,:,k+1)-u_ini(:,:,k-1))/(dzt(k+1)+dzt(k))
 !enddo
 !uz(:,:,1)=uz(:,:,2)
 !uz(:,:,nz)=uz(:,:,nz-1)
 !do k=1,nz
 !  aloc(1,js_pe:je_pe) = dyt(js_pe:je_pe)*uz(is_pe,js_pe:je_pe,k) *coriolis_t(is_pe,js_pe:je_pe)/grav/alpha*rho_0
 !  call pe0_recv_2D(nx,ny,aloc)
 !  cloc(1:ny,k) = aloc(1,1:ny)
 !enddo
 !bloc(1,:)=0.
 !do j=1,ny-1
 !   bloc(j+1,:)=bloc(j,:) + cloc(j,:)
 !   call pe0_bcast(bloc(j+1,:),nz)
 !enddo

 !open(27,file='namelist.in', status='OLD')
 !read(27,nml=nml_seed)
 !close(27)
 !call srand(seed)
 !print*,'seed = ',seed

 !do j=js_pe,je_pe
 !  do k=1,nz
 !   do i=is_pe,ie_pe
 !     temp(i,j,k,:)  = temp(i,j,k,:) + bloc(j,k) !+ 1e-3*rand()
 !   enddo
 !  enddo
 !enddo



 if (enable_rossmix) then
  do k=2,nphi-1
    do j=js_pe,je_pe
     do i=is_pe,ie_pe
      E_s(i,j,k,:,1:nmodes)= 1e-9*rand()
      E_l(i,j,k,:,1:nmodes)= 1e-9*rand()
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
 use main_module
 use config_module
 implicit none
 integer :: i,k
 real*8 :: zmean(js_pe:je_pe)
 if (enable_momentum_sources) then
    do k=1,nz
       zmean = 0
       do i=is_pe,ie_pe
        zmean  = zmean + u(i,js_pe:je_pe,k,tau)
       enddo
       call zonal_sum_vec2(zmean,je_pe-js_pe+1)
       do i=is_pe,ie_pe
        u_source(i,js_pe:je_pe,k)=-u_rest*zmean(js_pe:je_pe)/nx
       enddo
    enddo
 endif
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





