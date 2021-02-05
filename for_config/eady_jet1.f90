
!=======================================================================
!  channel with zonal mean restoring and linear drag
!=======================================================================

module config_module
 ! use this module only locally in this file
 implicit none
 !integer,parameter :: H_RES = 120    ! lateral grid points
 !integer,parameter :: V_RES = 40     ! vertical grid points
 integer,parameter :: X_RES = 60    ! lateral grid points
 integer,parameter :: Y_RES = 60*2  ! lateral grid points
 integer,parameter :: V_RES = 20     ! vertical grid points
 !
 !real*8,parameter :: Ri  = 20    ! Richardson number 
 !real*8,parameter :: CFL = 0.05  ! CFL number  
 !real*8,parameter :: delta = 0.02! aspect ratio
 !real*8,parameter :: Ek  = 0.01  ! Ekman number
 ! or
 real*8,parameter :: Ri  = 400    ! Richardson number 
 real*8,parameter :: CFL = 0.02   ! CFL number  
 real*8,parameter :: delta = 0.01 ! aspect ratio
 real*8,parameter :: Ek  = 0.050  ! Ekman number
 ! or
 !real*8,parameter :: Ri  = 5      ! Richardson number 
 !real*8,parameter :: CFL = 0.02   ! CFL number  
 !real*8,parameter :: delta = 0.02 ! aspect ratio
 !real*8,parameter :: Ek  = 0.01   ! Ekman number
 !
 real*8,parameter :: Ro  = sqrt( 1./Ri) ! Rossby number 
 real*8,parameter :: f0   = 1e-4        ! Coriolis freq.
 real*8,parameter :: H0   = 200.0       ! total depth
 real*8,parameter :: N0=f0/delta        ! stability freq.
 real*8,parameter :: U0 = Ro*N0*H0      ! mean flow
 real*8,parameter :: M0 = sqrt(f0*U0/H0)    ! meridional buoyancy gradient   
 real*8,parameter :: Lr = N0*H0/f0          ! Rossby radius
 real*8,parameter :: om_max = sqrt(5./54)/sqrt(1+Ri)*f0   ! growth rate of fastest growing mode
 real*8,parameter :: kmax   = 1./ ( sqrt(2./5.)*sqrt(1+Ri)*M0**2/f0**2*H0 )  ! wave number of fastest growing mode
 real*8,parameter :: Lx = 1*2/kmax  *3.14159265358979323846264338327950588
 real*8,parameter :: Ly = 4*2/kmax  *3.14159265358979323846264338327950588
 real*8,parameter :: jet_scale = 0.05*Ly
end module config_module



subroutine set_parameter
 ! ----------------------------------
 !       set here main parameter
 ! ----------------------------------
 use main_module
 use config_module
 use diagnostics_module
 use qg_module !!!!! new !!!!
 implicit none

 nx=X_RES; ny=Y_RES; nz = V_RES
     
 congr_epsilon = 1e-12
 congr_max_iterations = 5000
 enable_free_surface = .true. 
     
 enable_hydrostatic      = .true.
 enable_cyclic_x         = .true.
     
 !enable_superbee_advection     = .true.
 enable_upwind3_advection     = .true. !!!! new

 enable_explicit_vert_friction = .true.
 enable_biharmonic_friction    =  .true.

 enable_conserve_energy =  .false.
 coord_degree           =  .false.
 eq_of_state_type       =  1

 !!! new !!!
 enable_diag_qg_filter = .true.
 enable_qg_initialize = .true.
 !enable_qg_diag_Nsqr = .true.
 enable_qg_diag_ageo = .true.
 qg_Nsqr = N0**2
 qg_f0 = f0
 !!!!!!!!!!!
end subroutine set_parameter


subroutine set_grid
 use main_module   
 use config_module   
 use diagnostics_module
 use qg_module
 implicit none


 dzt = H0/nz
 dxt = Lx/nx
 dyt = Ly/ny

 dt_mom    = CFL/U0*dxt(is_pe)   !CFL=U*dt/dx
 dt_tracer = CFL/U0*dxt(is_pe)
 congr_epsilon = 1e-12 *(dxt(is_pe)/20e3)**2
 !a_h      = Ek*f0*dxt[0]**2 
 a_hbi    = Ek*f0*dxt(is_pe)**4 
 kappam_0 = Ek*f0*dzt(1)**2
 !M.k_h = Ek*f0*M.dxt[0]**2 
 !M.k_v = Ek*f0*M.dzt[0]**2 

 if (my_pe == 0) then
       print*,
       print*,' Lr = ',(Lr/1e3),' km'
       print*,' Ro = ',Ro
       print*,' Ri = ',Ri
       print*,' delta = ',delta
       print*,' Max. growth rate =' , om_max*86400,' 1/d'
       print*,' k_max Lr = ',kmax*Lr
       print*,' dx=',dxt(is_pe)/1e3,' km, dt= ',dt_mom,' s '
       print*,' CFL = ',U0*dt_mom/dxt(is_pe)
       print*,' CFL = ',om_max/kmax*dt_mom/dxt(is_pe)
       print*,' A_h = ',a_h,' m^2/s  Ek = ',Ek
       print*,' A_v = ',kappam_0,' m^2/s  Ek = ',Ek
 endif

 runlen =  86400.0*45
 enable_diag_ts_monitor = .true.; ts_monint = dt_mom!86400*0.25
 enable_diag_snapshots  = .true.; snapint   = 86400!dt_mom*0.25

 qg_snap_int = snapint
 qg_congr_epsilon = .001
 qg_congr_max_itts = 5000

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
 implicit none
 integer :: k,j,i
 real*8 :: alpha,get_drhodT
 real*8 :: u_ini(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) 
 real*8 :: aloc(nx,ny),bloc(ny,nz),cloc(ny,nz)
 real*8 :: uz(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz)
 real*8 :: ampA,ampB,fxa,d
 complex*16 :: c1,A,phiz

 alpha = get_drhodT(35d0,5d0,0d0) 
 salt = 35.0
 do k=1,nz
   temp(:,:,k,tau)=-N0**2*zt(k)/grav/alpha*rho_0*maskT(:,:,k)
 enddo

 u_ini = 0
 do j=js_pe,je_pe
  do k=1,nz
   u_ini(:,j,k)     = U0*exp(-(yt(j)-Ly/2)**2/jet_scale**2 )*cos(pi*zt(k)/H0)   !0.5*(1+tanh( (zt(k)+H0/2)/(0.2*H0) ))
  enddo
 enddo
 call border_exchg_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,u_ini)
 call setcyclic_xyz   (is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,u_ini)
 do k=2,nz-1
   uz(:,:,k) = (u_ini(:,:,k+1)-u_ini(:,:,k-1))/(dzt(k+1)+dzt(k))
 enddo
 uz(:,:,1)=uz(:,:,2)
 uz(:,:,nz)=uz(:,:,nz-1)


 u(:,:,:,tau)   = u_ini; 
 u(:,:,:,taum1) = u_ini; 
 u(:,:,:,taup1) = u_ini; 

 ! rho0 fu = -p_y, p_z = -g rho, rho0 f u_z = -g rho_y,  rho_y = - rho0 f u_z/g
 ! = alpha T_y
 do k=1,nz
   aloc(1,js_pe:je_pe) = dyt(js_pe:je_pe)*uz(is_pe,js_pe:je_pe,k)*f0/grav/alpha*rho_0
   call pe0_recv_2D(nx,ny,aloc)
   cloc(1:ny,k) = aloc(1,1:ny)
 enddo
 bloc(1,:)=0.
 do j=1,ny-1
    bloc(j+1,:)=bloc(j,:) + cloc(j,:)
    call pe0_bcast(bloc(j+1,:),nz)
 enddo

 
 call srand(1234)
 ampA = rand()
 ampB = rand()
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
    !fxa = sin(xt(i)*kmax)
    fxa = fxa*exp(-(yt(j)-Ly/2)**2/jet_scale**2 ) 
    fxa = fxa*abs(phiz)*rho_0/grav/alpha
    temp(i,j,k,:)  = temp(i,j,k,:) + bloc(j,k) + 0.2*fxa
   enddo
  enddo
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



