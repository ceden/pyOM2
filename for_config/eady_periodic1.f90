
!=======================================================================
!  double periodic domain with zonal mean restoring and linear drag
!=======================================================================

module config_module
 ! use this module only locally in this file
 implicit none
 integer,parameter :: H_RES = 60    ! lateral grid points
 integer,parameter :: V_RES = 40     ! vertical grid points
 !
 !real*8,parameter :: Ri  = 20    ! Richardson number 
 !real*8,parameter :: CFL = 0.02  ! CFL number  
 !real*8,parameter :: delta = 0.02! aspect ratio
 !real*8,parameter :: Ek  = 0.01  ! Ekman number
 ! or
 real*8,parameter :: Ri  = 400    ! Richardson number 
 real*8,parameter :: CFL = 0.01   ! CFL number  
 real*8,parameter :: delta = 0.01 ! aspect ratio
 real*8,parameter :: Ek  = 0.01  ! Ekman number
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
 real*8,parameter :: t_rest = om_max/5      ! restoring time scale
 real*8,parameter :: u_rest = om_max/5.   ! restoring time scale
 real*8 :: Lx,Ly
 real*8, allocatable :: U_0(:,:,:),T_0(:,:,:)
end module config_module



subroutine set_parameter
 ! ----------------------------------
 !       set here main parameter
 ! ----------------------------------
 use main_module
 use config_module
 use diagnostics_module
 use qg_module
 implicit none

 nx=H_RES; ny=H_RES; nz = V_RES
     
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
 coord_degree           =  .false.
 eq_of_state_type       =  1
 !enable_tempsalt_sources = .true.
 !enable_momentum_sources = .true.

 !!! new !!!
 !enable_diag_qg_filter = .true.
 !enable_qg_initialize = .true.
 !enable_qg_diag_Nsqr = .false.
 !qg_Nsqr = N0**2
 !qg_f0 = f0
 !qg_congr_epsilon = 1e-3
 !qg_press_epsilon = 1e-8
 !qg_congr_max_itts = 5000
 !!!!!!!!!!!

end subroutine set_parameter


subroutine set_grid
 use main_module   
 use config_module   
 use diagnostics_module
 use qg_module
 implicit none

 Lx = 4*2*pi/kmax   
 Ly = 4*2*pi/kmax   

 dzt = H0/nz
 dxt = Lx/nx
 dyt = Lx/ny

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
 enable_diag_snapshots  = .true.; snapint   = 86400*.1!dt_mom
 enable_diag_ts_monitor = .true.; ts_monint = snapint


 qg_snap_int = snapint
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
 use qg_module
 implicit none
 real*8 :: d,fxa
 real :: fxb
 integer :: k,j,i
 real*8 :: alpha,get_drhodT
 real*8 :: aloc(nx,ny),bloc(ny)
 real*4 :: r4

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

 fxa = M0**2*Ly/grav/alpha*rho_0
 do j=js_pe,je_pe
    T_0(:,j,:) = T_0(:,j,:) + fxa*sin(2*pi*yt(j)/Ly)
 enddo

 do k=1,nz
  temp(:,:,k,tau)   = T_0(:,:,k)
  temp(:,:,k,taum1) = T_0(:,:,k)
 enddo
 
 
 enable_qg_initialize = .true.
 enable_qg_diag_Nsqr = .false.
 enable_qg_diag_ageo = .false.  
 qg_congr_epsilon = 1e-4
 qg_press_epsilon = 1e-8
 qg_congr_max_itts = 5000
 qg_Nsqr = N0**2
 qg_f0 = f0
 call qg_filter_init
 enable_qg_initialize = .false.
 enable_diag_qg_filter = .false.
 do k=1,nz
  U_0(:,:,k) = u(:,:,k,tau)
 enddo



 ! perturbation buoyancy
 do k=1,nz
  do j=1,ny
   do i=1,nx
    call gauss3(r4); aloc(i,j)=r4
   enddo
  enddo
  call pe0_bcast(aloc,nx*ny)
  do j=js_pe,je_pe
    do i=is_pe,ie_pe
     !call random_number(fxb)
     !temp(:,j,k,tau) =T_0(:,j,k) + fxa*0.01*sin(kmax*xt)*sin(pi/Ly*yt(j))*sin(zt(k)*pi/H0)
     !temp(i,j,k,tau) =T_0(i,j,k) + fxa*(fxb-0.5)*0.01
     temp(i,j,k,tau) =T_0(i,j,k) + fxa*(aloc(i,j)-0.5)*0.01
    enddo    
  enddo    
 enddo    
 temp(:,:,:,taum1) = temp(:,:,:,tau)

end subroutine set_initial_conditions



subroutine GAUSS3(g1)
implicit none
      INTEGER gn,j
      REAL g1,r1
      g1=0.
      gn=100 ! increase gn to get better distribution
      DO j=1,gn
       call random_number(r1)
       g1=g1+r1
      ENDDO
      g1=(g1-.5*gn)/SQRT(FLOAT(gn))*SQRT(3.)*2.
end subroutine





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
        temp_source(i,js_pe:je_pe,k)=t_rest*(T_0(i,js_pe:je_pe,k)-zmean(js_pe:je_pe)/nx)
       enddo
    enddo
 endif
 if (enable_momentum_sources) then
    do k=1,nz
       zmean = 0
       do i=is_pe,ie_pe
        zmean  = zmean + u(i,js_pe:je_pe,k,tau)
       enddo
       call zonal_sum_vec2(zmean,je_pe-js_pe+1)
       do i=is_pe,ie_pe
        u_source(i,js_pe:je_pe,k)=u_rest*(U_0(i,js_pe:je_pe,k)-zmean(js_pe:je_pe)/nx)
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



