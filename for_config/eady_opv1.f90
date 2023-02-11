
!=======================================================================
!  channel with zonal mean restoring and linear drag
!=======================================================================

module config_module
 ! use this module only locally in this file
 implicit none
 integer,parameter :: fac=1
 integer,parameter :: X_RES = 30*fac  ! lateral grid points
 integer,parameter :: Y_RES = 30*fac  ! lateral grid points
 integer,parameter :: V_RES = 10*fac  ! vertical grid points
 
 real*8,parameter :: Ri  = 2000    ! Richardson number 
 real*8,parameter :: delta = 0.02 ! aspect ratio

 real*8,parameter :: Ro  = sqrt( 1./Ri) ! Rossby number 
 real*8,parameter :: f0   = 1e-4        ! Coriolis freq.
 real*8,parameter :: H0   = 1000.0       ! total depth
 real*8,parameter :: N0=f0/delta        ! stability freq.
 real*8,parameter :: U0 = Ro*N0*H0      ! mean flow
 real*8,parameter :: M0 = sqrt(f0*U0/H0)    ! meridional buoyancy gradient   
 real*8,parameter :: Lr = N0*H0/f0          ! Rossby radius
 real*8,parameter :: om_max = sqrt(5./54)/sqrt(1+Ri)*f0   ! growth rate of fastest growing mode
 real*8,parameter :: kmax   = 1./ ( sqrt(2./5.)*sqrt(1+Ri)*M0**2/f0**2*H0 )  ! wave number of fastest growing mode
 real*8,parameter :: Lx = 4*2/kmax  *3.14159265358979323846264338327950588
 real*8,parameter :: Ly = 4*2/kmax  *3.14159265358979323846264338327950588
 real*8,parameter :: jet_scale = 0.1*Ly
end module config_module



subroutine set_parameter
 ! ----------------------------------
 !       set here main parameter
 ! ----------------------------------
 use main_module
 use config_module
 use diagnostics_module
 use diag_opt_balance_module
 implicit none
 nx=X_RES; ny=Y_RES; nz = V_RES
 congr_max_iterations = 5000
 !enable_free_surface = .true. 
 congr_epsilon = 1e-9
 enable_hydrostatic      = .true.
 enable_cyclic_x         = .true.
 enable_conserve_energy =  .false.
 coord_degree           =  .false.
 eq_of_state_type       =  1
 
 AB_eps = 0.01
 dt_mom    = 1600./fac
 dt_tracer = dt_mom 
 runlen =  86400.0*4500
 enable_diag_ts_monitor = .true.; ts_monint = 86400
 enable_diag_snapshots  = .true.; snapint   = 86400*5.
 enable_diag_opt_balance = .true.; opt_balance_int = 86400*5.
 opt_balance_max_Itts = 5
 opt_balance_period  = 86400.*7.5
 opt_balance_average = 86400.*15
 opt_balance_temp_only = .true.
end subroutine set_parameter


subroutine set_grid
 use main_module   
 use config_module   
 use diagnostics_module
 implicit none

 dzt = H0/nz
 dxt = Lx/nx
 dyt = Ly/ny

 ! f U = A_h U/L^4  
 !a_hbi    = 0.05*f0*max(dxt(is_pe),dyt(js_pe))**4 
 !kappam_0 = 0.050*f0*dzt(1)**2
 !enable_superbee_advection     = .true.

 if (my_pe == 0) then
       print*,' '
       print*,' Lr = ',(Lr/1e3),' km'
       print*,' Ro = ',Ro
       print*,' Ri = ',Ri
       print*,' delta = ',delta
       print*,' Max. growth rate =' , om_max*86400,' 1/d'
       print*,' k_max Lr = ',kmax*Lr
       print*,' dx=',dxt(is_pe)/1e3,' km, dt= ',dt_mom,' s '
 endif

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
 implicit none
 integer :: k,j,i,n
 real*8 :: alpha,get_drhodT
 real*8 :: u_ini(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) 
 real*8 :: aloc(nx,ny),bloc(ny,nz),cloc(ny,nz)
 real*8 :: uz(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz)
 real*8 :: ampA,ampB,fxa,d
 complex*16 :: c1,A,phiz

 
 salt = 35.0
 do k=1,nz
   temp(:,:,k,tau) = N0**2*zt(k)*maskT(:,:,k)
 enddo

 u_ini = 0
 do j=js_pe,je_pe
  do k=1,nz
   u_ini(:,j,k)     = U0*exp(-(yt(j)-Ly/2)**2/jet_scale**2 )*cos(pi*zt(k)/H0)   
  enddo
 enddo
 call border_exchg_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,u_ini)
 call setcyclic_xyz   (is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,u_ini)
 do k=2,nz-1
   uz(:,:,k) = (u_ini(:,:,k+1)-u_ini(:,:,k-1))/(dzt(k+1)+dzt(k))
 enddo
 uz(:,:,1)=uz(:,:,2)
 uz(:,:,nz)=uz(:,:,nz-1)
 u(:,:,:,tau)   = u_ini; u(:,:,:,taum1) = u_ini; u(:,:,:,taup1) = u_ini; 

 
 do k=1,nz
   aloc(1,js_pe:je_pe) = dyt(js_pe:je_pe)*uz(is_pe,js_pe:je_pe,k)*f0
   call pe0_recv_2D(nx,ny,aloc)
   cloc(1:ny,k) = aloc(1,1:ny)
 enddo
 bloc(1,:)=0.
 do j=1,ny-1
    bloc(j+1,:)=bloc(j,:) + cloc(j,:)
    call pe0_bcast(bloc(j+1,:),nz)
 enddo

 
 !call srand(1234)
 ampA = 0!rand()
 ampB = 0.3!rand()
 if (my_pe==0) then
   print*,' AmpA = ',AmpA,' AmpB =',AmpB
 endif


 d=f0/N0/kmax
 fxa=(exp(H0/d)+exp(-H0/d))/(exp(H0/d)-exp(-H0/d))
 c1 = 1+0.25*(H0/d)**2-H0/d*fxa 
 c1=(sqrt( c1 )*d/H0+0.5)*U0
 A=(U0-c1)/U0*H0/d

 do i=is_pe,ie_pe
  do j=js_pe,je_pe
   do k=1,nz
    phiz=A/d*sinh(zt(k)/d)+cosh(zt(k)/d)/d
    fxa = ampA*sin(xt(i)*kmax)+ampB*cos(xt(i)*kmax)
    fxa = fxa*exp(-(yt(j)-Ly/2)**2/jet_scale**2 ) 
    fxa = fxa*abs(phiz)
    temp(i,j,k,:)  = temp(i,j,k,:) - bloc(j,k) - 0.2*fxa
   enddo
  enddo
 enddo
 
 alpha = get_drhodT(35d0,5d0,0d0) 
 temp = -temp*rho_0/grav/alpha
 
!if (.false.) then 
 call calc_initial_conditions

 
 call diag_opt_balance
 call diag_opt_write 
!call halt_stop(' ')


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
     call border_exchg_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,salt(:,:,:,n)) 
     call setcyclic_xyz   (is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,salt(:,:,:,n))
     call border_exchg_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,u(:,:,:,n)) 
     call setcyclic_xyz   (is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,u(:,:,:,n))
     call border_exchg_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,v(:,:,:,n)) 
     call setcyclic_xyz   (is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,v(:,:,:,n))
     call border_exchg_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,w(:,:,:,n)) 
     call setcyclic_xyz   (is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,v(:,:,:,n))
     call border_exchg_xy(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,psi(:,:,n)) 
     call setcyclic_xy   (is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,psi(:,:,n))

  enddo
!endif

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


