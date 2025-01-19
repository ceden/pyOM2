
! Atmosphere as ideal gas
! k runs here from p=0 to p=p_max, so top-down.

! initial conditions for temp as in Fig. 1 of Snyder et al 91, JAS, A comparison of 

module config_module
 ! use this module only locally in this file
 implicit none
 real*8,parameter :: fac=1.0
 real*8 :: f0=1e-4, Ly = 8000e3, Lx = 8000e3, Lz = 1e5 ! 1bar = 10^5 Pa
 real*8, parameter :: R_gas= 287.0 , cp=1e3, p_c= 1e5, u0 = 40.0, kappa = R_gas/cp
 real*8, parameter :: p_T = 1.5e4, dp_T = 0!*p_T
 real*8, parameter :: N_strato = 0.02, N_tropo = 0.01
end module config_module



subroutine set_parameter
 use main_module 
 use config_module
 use diagnostics_module
 use diag_opt_balance_module
 implicit none

 nx=int(50*fac); ny=int(50*fac); nz = int(20*fac)
 
 dt_mom    = 96./fac 
 dt_tracer = dt_mom 
 runlen =  86400.0*1000
 
 congr_max_iterations = 5000

 congr_epsilon = 1e-9
 enable_hydrostatic      = .true.
 enable_cyclic_x         = .true.
 enable_conserve_energy  = .false.
 coord_degree            = .false.
 eq_of_state_type        =  100
 
 enable_set_zero_at_surface = .false.
 AB_eps = 0.01

 enable_biharmonic_friction = .true.
 A_hbi = 1e11*(Lx/nx/40e3)**4   !1e12/fac**4 !   u_t = u/dx**4 , T = (dx/fac)**4
 enable_biharmonic_mixing = .true.
 K_hbi = A_hbi ! 1e12/fac**4
 !enable_superbee_advection     = .true.
 
 if (my_pe==0) then
  print*,' A_hbi = ',A_hbi,' m^4/s'
  print*,' K_hbi = ',K_hbi,' m^4/s'
 endif
 
 enable_diag_ts_monitor = .true.; ts_monint = dt_mom
 enable_diag_snapshots  = .true.; snapint   = 86400*2
  
 enable_diag_opt_balance = .true.; 
 opt_balance_int = snapint
 opt_balance_period  = 0.5* 2*pi/f0
 opt_balance_max_Itts = 1
 opt_balance_average_times = 5
 opt_balance_average = 2*pi*0.7/1e-4
 opt_balance_temp_only = .true. 

end subroutine set_parameter 


subroutine set_grid
 use main_module   
 use config_module   
 implicit none
 dxt=Lx/nx;dyt=Ly/ny;dzt=Lz/nz
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
 use diagnostics_module
 use diag_opt_balance_module
 implicit none
 integer :: i,j,k,n
 real*8 :: uz(nx,ny),Ty(nx,ny), dy, fxa, fxb,exner
 real*8,dimension(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) :: theta
 
 character*80 :: filename
 logical :: file_exists
 
 write(filename,'(a,i5,a)')  'restart_PE_',my_pe,'.dta'
 call replace_space_zero(filename)
 inquire ( FILE=filename, EXIST=file_exists )
 if (file_exists) return

 dy = dyt(js_pe)
 theta=0d0
 
 do k=1,nz
  do j=js_pe,je_pe
   fxb = p_T + dp_T*tanh( (yt(j)-Ly/2.)/1000e3 )
   if (zt(k) <fxb) fxa = zt(k)/fxb
   if (zt(k)>=fxb) fxa = 1.-(zt(k)-fxb)/(p_C-fxb)
   u(:,j,k,tau) = u0*exp(-(yt(j)-Ly/2)**2/1000e3**2 )*fxa*maskU(:,j,k)
  enddo
 enddo 
 call border_exchg_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,u(:,:,:,tau))
 call setcyclic_xyz   (is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,u(:,:,:,tau)) 
 u(:,:,:,taum1) = u(:,:,:,tau)
 
 ! p_z = - g \rho
 ! N^2 = g/theta dtheta/dz = g/theta dtheta/dp dp/dz = g^2/(theta*alpha) dtheta/dp

 theta(:,:,nz) = 300
 do k=nz-1,1,-1 
   do j=js_pe,je_pe
   exner = cp*(zt(k+1)/p_c)**kappa ! Exner function 
   fxb = p_T + dp_T*tanh( (yt(j)-Ly/2.)/1000e3 )
   if (zt(k) <fxb) fxa = N_strato**2
   if (zt(k)>=fxb) fxa = N_tropo**2
   theta(:,j,k) = theta(:,j,k+1) +  fxa*kappa*exner*theta(:,j,k+1)**2/zt(k+1)/grav**2*dzt(k+1)
  enddo 
 enddo
 theta = theta*maskT

 do k=1,nz
  
  if (k==1) then
   uz(is_pe:ie_pe,js_pe:je_pe) = (u(is_pe:ie_pe,js_pe:je_pe,k+1,tau)-u(is_pe:ie_pe,js_pe:je_pe,k  ,tau))/dzw(k)
  elseif (k==nz) then
   uz(is_pe:ie_pe,js_pe:je_pe) = (u(is_pe:ie_pe,js_pe:je_pe,k  ,tau)-u(is_pe:ie_pe,js_pe:je_pe,k-1,tau))/dzw(k-1)
  else
   uz(is_pe:ie_pe,js_pe:je_pe) = (u(is_pe:ie_pe,js_pe:je_pe,k+1,tau)-u(is_pe:ie_pe,js_pe:je_pe,k-1,tau))/(dzw(k)+dzw(k-1))
  endif
  call pe0_recv_2D(nx,ny,uz)
  call pe0_bcast(uz,nx*ny)
    
  Ty = 0d0
  fxa = f0*zt(k)**(1.-kappa)/(kappa*cp/p_c**kappa)  
  do j=2,ny
   Ty(:,j) = Ty(:,j-1) + uz(:,j)*dy*fxa
  enddo
      
  do j=js_pe,je_pe
   do i=is_pe,ie_pe
     fxa = 1*cos(pi*zt(k)/p_c)
     !fxa = 1*(1-zt(k)/p_c)    
     theta(i,j,k) = theta(i,j,k) + (Ty(i,j)  + 0.05*fxa*sin(2*2*pi*xt(i)/Lx))*maskT(i,j,k)
   enddo
  enddo  
 enddo
 
 call border_exchg_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,theta)
 call setcyclic_xyz   (is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,theta) 
 
 temp(:,:,:,tau)   = theta
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
    ! boundary exchange
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
 kbot=1
end subroutine set_topography


subroutine set_diagnostics
end subroutine set_diagnostics


subroutine set_particles
end subroutine set_particles
