
!=======================================================================
!  channel with zonal mean restoring and linear drag
!=======================================================================

module config_module
 ! use this module only locally in this file
 implicit none
 integer :: gridpoints =  128*2
 real*8,parameter :: Ld = 96e3, H0 = 4e3, f0=1e-4, beta_loc = 1.56e-11, Ustar = 0.6
 real*8,parameter :: Ek  = 0.050  , alpha = 1., Ljet = 600e3*2
 real*8 :: Lx,Ly,dx, T_rest = 1./(3*86400.)
 real*8, allocatable :: U0(:,:,:),T0(:,:,:)
end module config_module



subroutine set_parameter
 ! ----------------------------------
 !       set here main parameter
 ! ----------------------------------
 use main_module
 use config_module
 implicit none

 nx=gridpoints; ny=gridpoints*2; nz = gridpoints/4    
 congr_epsilon = 1e-12
 congr_max_iterations = 5000
 enable_free_surface = .true.      
 enable_hydrostatic      = .true.
 enable_cyclic_x         = .true.     
 enable_superbee_advection     = .true. 
 enable_explicit_vert_friction = .true.
 enable_biharmonic_friction    =  .true.
 enable_conserve_energy =  .false.
 coord_degree           =  .false.
 eq_of_state_type       =  1
 enable_tempsalt_sources = .true.
 enable_momentum_sources = .true. 
 !enable_bottom_friction = .true.
end subroutine set_parameter



subroutine set_grid
 use main_module   
 use config_module   
 use diagnostics_module
 implicit none
 real*8 :: S(nz),b0,b1,dz_loc,b(nz),z
 integer :: k 
 
 
 Lx = 8*Ld   
 Ly = Lx*2.  
 dx = Lx/nx
 dxt = dx
 dyt = dx

 ! vertical grid from dz = dx/S, S=N/f
 b1 = 0.5*(1+tanh(8*(-0.95)) )
 b0 = 1./(0.1-b1+0.5*(1+tanh(8*(1-0.95))))
 dz_loc = H0/nz
 do k=1,nz
  z = (k-1)*dz_loc + dz_loc/2. -H0
  b(k) = b0*(0.1*(1+z/H0 )-b1+0.5*(1+tanh(8*((1+z/H0 )-0.95))  ) )
 enddo
 do k=1,nz-1 
  S(k) = sqrt(  (b(k+1)-b(k))/dz_loc )
 enddo
 S(nz) = S(nz-1) 
 dzt = dx/S
 dzt = dzt/sum(dzt)*H0

 
 dt_mom    =  120.0*2 *(dx/3e3)
 dt_tracer = dt_mom
 AB_eps = 0.1
 congr_epsilon = 1e-6 *(dx/3e3)**2
 
 a_hbi    = Ek*f0*dx**4 
 kappam_0 = Ek*f0*( sum(dzt)/nz) **2
 r_bot = 1./(30.*86400.)*100./dzt(1) 
 
 runlen =  86400.0*365*5
 enable_diag_ts_monitor = .true.; ts_monint = dt_mom!86400*0.25
 !enable_diag_snapshots  = .true.; 
 enable_diag_parallel_snap=.true.
 snapint   =  86400*2.!dt_mom*0.25

end subroutine set_grid


subroutine set_coriolis
 use main_module   
 use config_module
 implicit none
 integer :: j
 do j=js_pe-onx,je_pe+onx
  coriolis_t(:,j) = f0 + beta_loc*(yt(j)-Ly/2)
 enddo
end subroutine set_coriolis


subroutine set_initial_conditions
 use main_module
 use config_module
 implicit none
 real*8 :: alpha2,b,zbar,b0,b1,get_drhodT
 real*8 :: aloc(nx,ny),p0(ny,nz),bs(ny,nz),r(nx)
 integer :: i,j,k
 
 allocate( U0(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) )
 allocate( T0(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) )
 U0 = 0; T0 = 0
    
 alpha2 = get_drhodT(35d0,5d0,0d0)  
 ! b = - g/rho0 alpha2 T ,  T = - rho_0/(g alpha2) b 
 
 ! z=0 -> zbar=1,  b= (f0*Ld)**2/H0*b0*(0.1-b1+0.5*(1+tanh(8*(1-0.95))  ) )
 ! z=-h-> bar=0 , b=(f0*Ld)**2/H0*b0*(-b1+0.5*(1+tanh(8*(-0.95))  ) )
 b1 = 0.5*(1+tanh(8*(-0.95)) )
 b0 = 1./(0.1-b1+0.5*(1+tanh(8*(1-0.95))))
 
 ! vertical shape of background buoyancy and zonal mean zonal flow
 do k=1,nz
   zbar = 1+zt(k)/H0  
   b = (f0*Ld)**2/H0*b0*(0.1*zbar-b1+0.5*(1+tanh(8*(zbar-0.95))  ) )
   T0(:,:,k) = -rho_0/grav/alpha2*b*maskT(:,:,k)
   U0(:,:,k) =  Ustar*beta_loc*Ld**2 *(  2*(alpha-1)*(b*H0/(f0*Ld)**2)**3+(3-2*alpha)*(b*H0/(f0*Ld)**2)**2 )
 enddo 
 
 ! meridional shape of zonal mean zonal flow
 do j=js_pe-onx,je_pe+onx
  if ( abs(yu(j)-Ly/2.)/Ljet < 0.5) then 
    U0(:,j,:) = U0(:,j,:)*cos( (yu(j)-Ly/2.)/Ljet*pi) 
  else 
    U0(:,j,:) = 0.
  endif
 enddo
 u(:,:,:,tau)     = U0
 u(:,:,:,taum1)   =  u(:,:,:,tau)
  
 ! pressure p0 from integrated discrete geostrophic balance
 ! fu = -p_y,  p(j+1) = p_j -  (f_j u_j + f_j+1 u_j+1)dx/2  
 do k=1,nz
  aloc(is_pe:ie_pe,js_pe:je_pe) = coriolis_t(is_pe:ie_pe,js_pe:je_pe)* &
                                          U0(is_pe:ie_pe,js_pe:je_pe,k)
  call pe0_recv_2D(nx,ny,aloc)
  p0(1,k)=0.
  do j=1,ny-1
    p0(j+1,k) = p0(j,k) - (aloc(1,j) + aloc(1,j+1) )*dx/2.  
  enddo
  call pe0_bcast(p0(:,k),ny)
 enddo

 ! perturbation buoyancy bs from discrete hydrostatic relation
 !p_z = b,   (p_k+1-p_k  )2/dz -b_k = b_k+1
 bs(:,1) = 0.
 do k=1,nz-1
  bs(:,k+1) = (p0(:,k+1)-p0(:,k))*2/dzw(k) - bs(:,k)
 enddo

 ! total background temperature
 do k=1,nz
  do j=js_pe,je_pe 
   T0(:,j,k) = T0(:,j,k) - rho_0/(grav*alpha2)*bs(j,k)*maskT(:,j,k)
  enddo
 enddo 
 call border_exchg_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,T0) 
 call setcyclic_xyz   (is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,T0)
 temp(:,:,:,tau)  = T0 
 temp(:,:,:,taum1) =  temp(:,:,:,tau)
  
 ! small perturbation in v  
 do k=1,nz
  do j=js_pe,je_pe
   do i=is_pe,ie_pe    
     v(i,j,k,tau) =  Ustar*beta_loc*Ld**2/100.*sin( xt(i)/Lx*12*pi)*cos(yu(j)/Ly*pi)*cos(zt(k)/H0*pi)      
   enddo
  enddo
 enddo 
 
 do k=1,nz 
  do j=js_pe,je_pe
   call random_number(r)
   call pe0_bcast(r,nx)
   v(is_pe:ie_pe,j,k,tau) = v(is_pe:ie_pe,j,k,tau) + Ustar*beta_loc*Ld**2/500*(r(is_pe:ie_pe)-0.5)
  enddo
 enddo
 
 call border_exchg_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,v(:,:,:,tau)) 
 call setcyclic_xyz   (is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,v(:,:,:,tau))
 v(:,:,:,taum1)   =  v(:,:,:,tau)
 
end subroutine set_initial_conditions



subroutine set_forcing
 use main_module
 use config_module
 implicit none
 integer :: i,k
 real*8 :: zmean(js_pe:je_pe)

 if (enable_tempsalt_sources) then
    do k=1,nz
       zmean = 0
       do i=is_pe,ie_pe
        zmean  = zmean + temp(i,js_pe:je_pe,k,tau)
       enddo
       call zonal_sum_vec2(zmean,je_pe-js_pe+1)
       do i=is_pe,ie_pe
        temp_source(i,js_pe:je_pe,k)=T_rest*(T0(i,js_pe:je_pe,k)-zmean(js_pe:je_pe)/nx)
       enddo
    enddo
    !call border_exchg_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,temp_source) 
    !call setcyclic_xyz   (is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,temp_source)
 endif
 if (enable_momentum_sources) then
    do k=1,nz
       zmean = 0
       do i=is_pe,ie_pe
        zmean  = zmean + u(i,js_pe:je_pe,k,tau)
       enddo
       call zonal_sum_vec2(zmean,je_pe-js_pe+1)
       do i=is_pe,ie_pe
        u_source(i,js_pe:je_pe,k)=T_rest*(U0(i,js_pe:je_pe,k)-zmean(js_pe:je_pe)/nx)
       enddo
    enddo
    !call border_exchg_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,u_source) 
    !call setcyclic_xyz   (is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,u_source)
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



