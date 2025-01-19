
! testing optimal balance with time averaging implementation
! two or one unstable jets, hydrostatic setup 

module config_module
 ! use this module only locally in this file
 implicit none
 real*8 :: f0,N0,Lx,Ly,Lz,fac=1,u0=0.2
end module config_module


subroutine set_parameter
 use main_module   
 use diagnostics_module
 use config_module
 use diag_opt_balance_module 
 implicit none
 
  !enable_cyclic_y      =.true.
  enable_cyclic_x      =.true.
  
  Lx = 500e3; Ly=Lx; Lz = 500.0
  f0 = 1e-4; N0 = 50*f0; dt_mom = 400/fac! 4./N0
  
  nx = int(100*fac); ny = int(100*fac); nz = int(10*fac)
  if (.not.enable_cyclic_y) then
   ny = ny/2
   Ly = Ly/2.
  endif 
  
  dt_tracer = dt_mom
  enable_conserve_energy = .false. 
  coord_degree           = .false. 
  enable_hydrostatic     = .true.
  eq_of_state_type = 1 
  congr_max_iterations = 5000
  congr_epsilon = 1e-12!1e-9
  !congr_epsilon_non_hydro=   1d-6
  AB_eps = 0.01  
  
  !enable_biharmonic_friction = .true.
  !A_hbi = 2e10/fac**4 
  
  enable_diag_snapshots  = .true.; snapint  = 15*86400.
  enable_diag_ts_monitor = .true.; ts_monint = dt_mom
  runlen = 40000*snapint
  
  enable_diag_opt_balance = .true.; opt_balance_int = snapint 
  opt_balance_period  = 2* 2*pi/f0
  opt_balance_max_Itts = 2
  
  opt_balance_average_times = 5
  opt_balance_average = 2*pi*0.7/1e-4
  opt_balance_temp_only = .true. 
  
end subroutine set_parameter


subroutine set_grid
 use main_module
 use config_module
 implicit none
 dxt(:)=Lx/nx
 dyt(:)=Ly/ny
 dzt(:)=Lz/nz
end subroutine set_grid


subroutine set_coriolis
 use main_module   
 use config_module
 implicit none
 coriolis_t = f0
end subroutine set_coriolis


subroutine set_initial_conditions
 use main_module  
 use config_module 
 use diag_opt_balance_module 
 use diagnostics_module 
 implicit none
 integer :: i,j,k,n
 real*8 :: alpha,get_drhodT
 real*8 :: bs(nx,ny)
 real*8 :: x(nx),y(ny),z(nz)
 
 do i=1,nx
   x(i)=(i-1)*dxt(is_pe)
 enddo 
 do i=1,ny
   y(i)=(i-1)*dyt(js_pe)
 enddo 
 do i=1,nz
   z(i)=(i-1)*dzt(1)
 enddo 
 
 alpha = get_drhodT(35d0,5d0,0d0) 

 if (enable_cyclic_y) then
  do k=1,nz
   do j=js_pe,je_pe
    do i=is_pe,ie_pe  
     u(i,j,k,tau) = u0*( ( exp( -(y(j)-3*Ly/4.)**2/(Lx*0.02)**2 ) &
                          -exp( -(y(j)-1*Ly/4.)**2/(Lx*0.02)**2 )  )*cos(z(k)/Lz*pi) &
          + 0.05*sin(x(i)/Lx*10*pi)*sin(y(j)/Ly*2*pi)*cos(z(k)/Lz*pi) )     
    enddo
   enddo
  enddo
 else
  do k=1,nz
   do j=js_pe,je_pe
    do i=is_pe,ie_pe  
     u(i,j,k,tau) = u0*(  -exp( -(y(j)-Ly/2.)**2/(Lx*0.02)**2 )*cos(z(k)/Lz*pi) &
                     + 0.05*sin(x(i)/Lx*10*pi)*sin(y(j)/Ly*2*pi)*cos(z(k)/Lz*pi) )     
    enddo
   enddo
  enddo
 endif
 call border_exchg_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,u(:,:,:,tau)) 
 call setcyclic_xyz   (is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,u(:,:,:,tau))
 u(:,:,:,taum1) = u(:,:,:,tau)

 
 do k=1,nz
  if (k==1) then
   bs(is_pe:ie_pe,js_pe:je_pe) = (u(is_pe:ie_pe,js_pe:je_pe,k+1,tau)-u(is_pe:ie_pe,js_pe:je_pe,k  ,tau))/(dzt(1))
  elseif (k==nz) then
   bs(is_pe:ie_pe,js_pe:je_pe) = (u(is_pe:ie_pe,js_pe:je_pe,k  ,tau)-u(is_pe:ie_pe,js_pe:je_pe,k-1,tau))/(dzt(1))
  else
   bs(is_pe:ie_pe,js_pe:je_pe) = (u(is_pe:ie_pe,js_pe:je_pe,k+1,tau)-u(is_pe:ie_pe,js_pe:je_pe,k-1,tau))/(2*dzt(1))
  endif
  bs = -dyt(js_pe)*bs*f0 
  call pe0_recv_2D(nx,ny,bs)
  call pe0_bcast(bs,nx*ny)
  
  ! 0 = -p_z + b, 0 = - p_y - fu,  b_y = -fu_z,   B = N0^2*z
  do j=2,ny
    bs(:,j) = bs(:,j) + bs(:,j-1)
  enddo
  ! rho = alpha T , b = - g/rho0 rho, temp  = - b*rho0/g/alpha,  b = - T g/rho0 alpha 
  do j=js_pe,je_pe
   do i=is_pe,ie_pe
    temp(i,j,k,tau) = -(N0**2*zt(k) + bs(i,j)) *rho_0/grav/alpha*maskT(i,j,k)
   enddo
  enddo 
 enddo
 call border_exchg_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,temp(:,:,:,tau))
 call setcyclic_xyz   (is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,temp(:,:,:,tau))
 temp(:,:,:,taum1) = temp(:,:,:,tau) 
 salt = 35d0
 
 if (enable_diag_opt_balance.and. .true.) then
 
  call calc_initial_conditions
  call diag_opt_balance
  call diag_opt_write 


  do i=is_pe,ie_pe
   do j=js_pe,je_pe
    do k=1,nz
     u(i,j,k,:) = u_bal(i,j,k)
     v(i,j,k,:) = v_bal(i,j,k)
     w(i,j,k,:) = w_bal(i,j,k)
     temp(i,j,k,:)  = temp_bal(i,j,k) + temp_ave(i,j,k)
     psi(i,j,:)  = psi_bal(i,j)
    enddo
   enddo
  enddo

  do n=1,3
     call border_exchg_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,temp(:,:,:,n)) 
     call setcyclic_xyz   (is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,temp(:,:,:,n))
     call border_exchg_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,u(:,:,:,n)) 
     call setcyclic_xyz   (is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,u(:,:,:,n))
     call border_exchg_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,v(:,:,:,n)) 
     call setcyclic_xyz   (is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,v(:,:,:,n))
     call border_exchg_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,w(:,:,:,n)) 
     call setcyclic_xyz   (is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,w(:,:,:,n))
     call border_exchg_xy(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,psi(:,:,n)) 
     call setcyclic_xy   (is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,psi(:,:,n))
  enddo



 endif
 
end subroutine set_initial_conditions


subroutine set_forcing
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