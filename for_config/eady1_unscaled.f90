
! two balanced jets with small perturbation get unstable
! testing optimal balance with time averaging

module config_module
 ! use this module only locally in this file
 implicit none
 integer,parameter :: fac=1
 real*8,parameter :: Ri  = 2000    ! Richardson number 
 real*8,parameter :: delta = 0.02 ! aspect ratio
 real*8,parameter :: H0   = 1000.0       ! total depth
 real*8 :: Lx,Ly,Lz,f0,N0,U0_loc,M0,kmax,jet_scale
end module config_module



subroutine set_parameter
 use main_module 
 use config_module
 use diagnostics_module
 use diag_opt_balance_module
 implicit none

 f0 = 1e-4
 N0 = f0/delta
 U0_loc = sqrt(1./Ri)*N0*H0   
 M0 = sqrt(f0*U0_loc/H0)  
 kmax   = 1./ ( sqrt(2./5.)*sqrt(1+Ri)*M0**2/f0**2*H0 )  
 Lx = 4*2/kmax  *pi
 Ly = 4*2/kmax  *pi * 2
 Lz = H0
 jet_scale = 0.1*Ly/2
 


 nx=30*fac; ny=30*2*fac; nz = 10*fac
 congr_max_iterations = 5000
 enable_free_surface = .true. 
 congr_epsilon = 1e-9
 enable_hydrostatic      = .true.
 enable_cyclic_x         = .true.
 enable_cyclic_y         = .true.
 enable_conserve_energy =  .false.
 coord_degree           =  .false.
 eq_of_state_type       =  1
 
 AB_eps = 0.01
 dt_mom    = 200./fac
 dt_tracer = dt_mom 
 runlen =  86400.0*4500
 enable_diag_ts_monitor = .true.; ts_monint = 86400*5
 enable_diag_snapshots  = .true.; snapint   = 86400*5.
 
 enable_diag_opt_balance = .true.; opt_balance_int = 86400*5.
 opt_balance_max_Itts = 3
 opt_balance_period  = 10./f0
 opt_balance_average = 2*86400.
 opt_balance_temp_only = .true. 
 opt_balance_average_times = 5
 

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
 integer :: i,j,k
 real*8 :: aloc(nx,ny),bloc(nx,ny)
 real*8 :: alpha,get_drhodT
  
 do k=1,nz
   do j=js_pe,je_pe
    do i=is_pe,ie_pe  
     u(i,j,k,tau) = U0_loc*( exp(-(yt(j)-1*Ly/4)**2/jet_scale**2) &
                        -exp(-(yt(j)-3*Ly/4)**2/jet_scale**2 ) &
              + 0.05*sin(xt(i)/Lx*4*pi)*sin(yt(j)/Ly*2*pi)  ) *cos(pi*zt(k)/H0)      
   enddo
  enddo
 enddo 
 call border_exchg_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,u(:,:,:,tau))
 call setcyclic_xyz   (is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,u(:,:,:,tau))

 do k=2,nz-1
   do j=js_pe,je_pe
    do i=is_pe,ie_pe  
     aloc(i,j) = -Ly/ny*(u(i,j,k+1,tau)-u(i,j,k-1,tau))/(2*Lz/nz) *f0
    enddo
   enddo     
   call pe0_recv_2D(nx,ny,aloc)  
   bloc(:,1)=0.
   do j=1,ny-1    
      bloc(:,j+1) = bloc(:,j) + aloc(:,j)
   enddo 
   call pe0_bcast(bloc,nx*ny)  
   temp(is_pe:ie_pe,js_pe:je_pe,k,tau)= bloc(is_pe:ie_pe,js_pe:je_pe)
 enddo  
 temp(:,:,1,tau) = temp(:,:,2,tau) 
 temp(:,:,nz,tau) = temp(:,:,nz-1,tau) 


 do k=1,nz
   temp(:,:,k,tau) = temp(:,:,k,tau)  + N0**2*zt(k)*maskT(:,:,k)
 enddo
 
 salt = 35.0
 alpha = get_drhodT(35d0,5d0,0d0) 
 temp = -temp*rho_0/grav/alpha
 
 
 if (enable_diag_opt_balance .and. .true. ) then
 
  call diag_opt_balance()
  
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
  
 
 
  do i=1,3
    ! boundary exchange
     call border_exchg_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,temp(:,:,:,i)) 
     call setcyclic_xyz   (is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,temp(:,:,:,i))
     call border_exchg_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,salt(:,:,:,i)) 
     call setcyclic_xyz   (is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,salt(:,:,:,i))
     call border_exchg_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,u(:,:,:,i)) 
     call setcyclic_xyz   (is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,u(:,:,:,i))
     call border_exchg_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,v(:,:,:,i)) 
     call setcyclic_xyz   (is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,v(:,:,:,i))
     call border_exchg_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,w(:,:,:,i)) 
     call setcyclic_xyz   (is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,v(:,:,:,i))
     call border_exchg_xy(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,psi(:,:,i)) 
     call setcyclic_xy   (is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,psi(:,:,i))

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
