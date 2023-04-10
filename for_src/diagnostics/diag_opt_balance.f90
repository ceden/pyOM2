!=======================================================================
! optimal balance with time averaging
! to diagnose fast and slow mode
!=======================================================================


module diag_opt_balance_module
 implicit none
 integer :: opt_balance_max_Itts = 5
 real*8 :: opt_balance_period  = 86400.*5
 real*8 :: opt_balance_average = 86400.*5
 integer :: opt_balance_average_times = 1
 real*8 :: opt_balance_tol     = 1d-9
 logical :: opt_balance_temp_only = .false.
 
 real*8 :: norm_diff
 real*8 :: opt_rho,opt_A_ab,opt_B_ab
 real*8,allocatable, dimension(:,:,:) :: u_base,v_base,w_base,temp_base,salt_base,rho_base
 real*8,allocatable, dimension(:,:)   :: psi_base,psi_bal
 real*8,allocatable, dimension(:,:,:) :: u_bal,v_bal,w_bal,temp_bal,salt_bal,rho_bal
 real*8,allocatable, dimension(:,:,:) :: temp_ave,salt_ave
 real*8,allocatable, dimension(:,:,:,:) :: dtemp_ave,dsalt_ave
end module diag_opt_balance_module



subroutine diag_opt_balance
!=======================================================================
! 
!=======================================================================
 use main_module
 use diag_opt_balance_module
 implicit none
 integer :: n,m,n_end_ramp,n_end_ave
 real*8,allocatable, dimension(:,:,:,:) :: u_bak,v_bak,w_bak,du_bak,dv_bak
 real*8,allocatable, dimension(:,:,:)   :: psi_bak,dpsi_bak
 real*8,allocatable, dimension(:,:,:,:) :: temp_bak,salt_bak,dtemp_bak,dsalt_bak
 real*8,allocatable, dimension(:,:,:,:) :: rho_bak
 logical, save :: first = .true.
 real*8,allocatable, dimension(:,:,:) :: u_loc,v_loc,w_loc,temp_loc,salt_loc,rho_loc
 real*8,allocatable, dimension(:,:)   :: psi_loc
 real*8,allocatable, dimension(:,:,:) :: u_save,v_save,w_save,temp_save,salt_save,rho_save
 real*8,allocatable, dimension(:,:)   :: psi_save
 integer :: o_tau,o_taum1,o_taup1
 real*8 :: fxa,fxb,fxc
 integer :: i,j,k
 
 if (my_pe==0) print*,'Entering optimal balance procedure '
 
 if (first) then
 
  if (my_pe==0) print*,' opt_balance_max_Itts = ',opt_balance_max_Itts 
  if (my_pe==0) print*,' opt_balance_period  = ',opt_balance_period, ' s'
  if (my_pe==0) print*,' opt_balance_average = ',opt_balance_average, ' s'
  if (my_pe==0) print*,' opt_balance_tol     = ',opt_balance_tol

  allocate( u_base(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); u_base = 0
  allocate( v_base(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); v_base = 0
  allocate( w_base(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); w_base = 0
  allocate( temp_base(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); temp_base = 0
  allocate( salt_base(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); salt_base = 0
  allocate( rho_base(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); rho_base = 0
  allocate( psi_base(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx) ); psi_base = 0

  allocate( u_bal(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); u_bal = 0
  allocate( v_bal(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); v_bal = 0
  allocate( w_bal(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); w_bal = 0
  allocate( temp_bal(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); temp_bal = 0
  allocate( salt_bal(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); salt_bal = 0
  allocate( rho_bal(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) );   rho_bal = 0
  allocate( psi_bal(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx)    );   psi_bal = 0
  
  allocate( temp_ave(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); temp_ave=0
  allocate( salt_ave(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); salt_ave=0
  allocate( dtemp_ave(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz,3) ); dtemp_ave=0
  allocate( dsalt_ave(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz,3) ); dsalt_ave=0
 
 endif
 
 
 ! store complete present state
 o_tau=tau; o_taum1=taum1; o_taup1=taup1
 allocate( u_bak(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz,3) ); u_bak = u
 allocate( v_bak(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz,3) ); v_bak = v
 allocate( w_bak(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz,3) ); w_bak = w 
 allocate( du_bak(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz,3) ); du_bak = du
 allocate( dv_bak(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz,3) ); dv_bak = dv
 
 allocate( psi_bak(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,3) ); psi_bak = psi
 allocate( dpsi_bak(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,3) ); dpsi_bak = dpsi

 allocate(temp_bak(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz,3) );  temp_bak = temp
 allocate(dtemp_bak(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz,3) ); dtemp_bak = dtemp
 allocate(salt_bak(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz,3) );  salt_bak = salt
 allocate(dsalt_bak(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz,3) ); dsalt_bak = dsalt
 allocate( rho_bak(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz,3) );  rho_bak = rho

 ! allocate local arrays for model state
 allocate( u_loc(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); u_loc=0
 allocate( v_loc(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); v_loc=0
 allocate( w_loc(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); w_loc=0
 allocate( temp_loc(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); temp_loc=0
 allocate( salt_loc(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); salt_loc=0
 allocate( rho_loc(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); rho_loc=0
 allocate( psi_loc(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx) ); psi_loc=0
 
 allocate( u_save(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); u_save=0
 allocate( v_save(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); v_save=0
 allocate( w_save(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); w_save=0
 allocate( temp_save(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); temp_save=0
 allocate( salt_save(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); salt_save=0
 allocate( rho_save(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); rho_save=0
 allocate( psi_save(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx) ); psi_save=0


 ! some checks
 if (enable_streamfunction)         call halt_stop('No streamfct. possible in opt_balance')
 if (enable_conserve_energy)        call halt_stop('No energy conservation in opt_balance')
 if (.not. enable_hydrostatic)      call halt_stop('No non-hydrostatic model in opt_balance')
 if (coord_degree)                  call halt_stop('No spherical coordinates in opt_balance')
 if (.not. enable_AB_time_stepping) call halt_stop('Need AB time stepping in opt balance')
 if (enable_tempsalt_sources)       call halt_stop('No source terms in opt_balance')
   
   
 ! calculate horizontally averaged temperature and salinity  
 do k=1,nz
  fxa = 0.; fxb=0; fxc=0
  do j=js_pe,je_pe
   do i=is_pe,ie_pe
     fxa = fxa + area_t(i,j)*maskT(i,j,k)
     fxb = fxb + temp(i,j,k,tau)*area_t(i,j)*maskT(i,j,k)
     fxc = fxc + salt(i,j,k,tau)*area_t(i,j)*maskT(i,j,k)
   enddo
  enddo
  call global_sum(fxa)
  call global_sum(fxb)
  call global_sum(fxc)
  temp_ave(:,:,k) = fxb/fxa
  salt_ave(:,:,k) = fxc/fxa
 enddo
 
 
 ! remove horizontally averaged temperature and salinity from model state
 temp(:,:,:,tau) = temp(:,:,:,tau) - temp_ave
 salt(:,:,:,tau) = salt(:,:,:,tau) - salt_ave

 ! save current model state 
 u_bal = u(:,:,:,tau); v_bal = v(:,:,:,tau); w_bal = w(:,:,:,tau)
 temp_bal = temp(:,:,:,tau); salt_bal = salt(:,:,:,tau); rho_bal = rho(:,:,:,tau)
 psi_bal = psi(:,:,tau)
 
 ! Calculate number of iteration steps for turning off/on the non linear terms
 n_end_ramp = int(opt_balance_period / dt_tracer) + 1
 n_end_ave  = int(opt_balance_average / dt_tracer) + 1

 ! set well defined model state
 do i=1,3
  u(:,:,:,i) = u_bal; v(:,:,:,i) = v_bal; w(:,:,:,i) = w_bal; psi(:,:,i) = psi_bal
  temp(:,:,:,i) = temp_bal; salt(:,:,:,i) = salt_bal; rho(:,:,:,i) = rho_bal
 enddo

 ! save base point by time averaging in linear model
 if (my_pe==0)  print '(a,i5,a,i5,a)',' integrating linear model for ',&
                       n_end_ave*opt_balance_average_times,' timesteps in ', &
                       opt_balance_average_times,' chunks'

 do m=1,opt_balance_average_times
  call diag_opt_time_ave(n_end_ave,u_base,v_base,w_base,temp_base,salt_base,rho_base,psi_base)
  do i=1,3
     u(:,:,:,i) = u_base;v(:,:,:,i) = v_base; w(:,:,:,i) = w_base; psi(:,:,i) = psi_base
     temp(:,:,:,i) = temp_base; salt(:,:,:,i) = salt_base; rho(:,:,:,i)  = rho_base
  enddo 
 enddo

 ! restore model state
 do i=1,3
  u(:,:,:,i) = u_bal; v(:,:,:,i) = v_bal; w(:,:,:,i) = w_bal; psi(:,:,i) = psi_bal
  temp(:,:,:,i) = temp_bal; salt(:,:,:,i) = salt_bal; rho(:,:,:,i) = rho_bal
 enddo
 
 norm_diff = 1d0

 ! Start optimal balance iteration
 outer: do n=1,opt_balance_max_Itts
 
    if (my_pe==0)  print '(a,i5)',' start optimal balance iteration:',n

    ! integrate backward to the linear end
    if (my_pe==0)  print '(a,i5,a)',' integrating backward to linear end for ',n_end_ramp,' timesteps '
    call diag_opt_backward_integration(n_end_ramp)

    ! project on geostrophic mode at linear end
    if (my_pe==0)  print '(a,i5,a,i5,a)',' integrating linear model for ',&
                       n_end_ave*opt_balance_average_times,' timesteps in ', &
                       opt_balance_average_times,' chunks'
    do m=1,opt_balance_average_times
     call diag_opt_time_ave(n_end_ave,u_loc,v_loc,w_loc,temp_loc,salt_loc,rho_loc,psi_loc)
     do i=1,3
      u(:,:,:,i) = u_loc; v(:,:,:,i) = v_loc; w(:,:,:,i) = w_loc; psi(:,:,i) = psi_loc
      temp(:,:,:,i) = temp_loc; salt(:,:,:,i) = salt_loc; rho(:,:,:,i) = rho_loc
     enddo
    enddo
    
    ! integrate forward to the non linear end
    if (my_pe==0)  print '(a,i5,a)',' integrating forward to non-linear end for ',n_end_ramp,' timesteps '
    call diag_opt_forward_integration(n_end_ramp)

    ! here we need to save the model state
    u_save = u(:,:,:,1); v_save = v(:,:,:,1); w_save = w(:,:,:,1); psi_save = psi(:,:,1)
    temp_save = temp(:,:,:,1); salt_save = salt(:,:,:,1); rho_save = rho(:,:,:,1)
    
    ! apply boundary condition at non linear end
    if (my_pe==0)  print '(a,i5,a,i5,a)',' integrating linear model for ',&
                       n_end_ave*opt_balance_average_times,' timesteps in ', &
                       opt_balance_average_times,' chunks'
    do m=1,opt_balance_average_times
     call diag_opt_time_ave(n_end_ave,u_loc,v_loc,w_loc,temp_loc,salt_loc,rho_loc,psi_loc)
     do i=1,3
      u(:,:,:,i) = u_loc; v(:,:,:,i) = v_loc; w(:,:,:,i) = w_loc; psi(:,:,i) = psi_loc
      temp(:,:,:,i) = temp_loc; salt(:,:,:,i) = salt_loc; rho(:,:,:,i) = rho_loc
     enddo
    enddo

    ! exchange base point
    do i=1,3
     u(:,:,:,i) = u_save - u_loc + u_base
     v(:,:,:,i) = v_save - v_loc + v_base
     w(:,:,:,i) = w_save- w_loc + w_base
     psi(:,:,i) = psi_save - psi_loc + psi_base
     temp(:,:,:,i) = temp_save - temp_loc + temp_base
     salt(:,:,:,i) = salt_save - salt_loc + salt_base
     rho(:,:,:,i)  = rho_save  - rho_loc  + rho_base
    enddo
    
    call diag_opt_balance_norm(u_bal,v_bal)

    ! update reference state
    u_bal = u(:,:,:,1); v_bal = v(:,:,:,1); w_bal = w(:,:,:,1); psi_bal = psi(:,:,1)
    temp_bal = temp(:,:,:,1); salt_bal = salt(:,:,:,1); rho_bal = rho(:,:,:,1)
    
    ! Check tolerance criterion
    if (n > 1) then
      if (my_pe==0)  print '(a,i5,a,ES15.3E3,a,ES15.3E3)', & 
              ' norm of difference to n= ',n - 1,' is ', norm_diff,' / ',opt_balance_tol
      if (norm_diff < opt_balance_tol) exit outer
    end if
 end do outer

 if (my_pe==0) print*,'Finishing optimal balance procedure '
 
 deallocate(u_save,v_save,w_save,temp_save,salt_save,rho_save,psi_save)
 deallocate(u_loc,v_loc,w_loc,temp_loc,salt_loc,rho_loc,psi_loc)
 
 ! restore complete present model state
 tau=o_tau; taum1=o_taum1; taup1=o_taup1
 dt_tracer = +abs(dt_tracer); dt_mom = + abs(dt_mom)
 u=u_bak;v=v_bak;w=w_bak;du=du_bak;dv=dv_bak;psi=psi_bak;dpsi=dpsi_bak
 temp=temp_bak;salt=salt_bak;dtemp=dtemp_bak;dsalt=dsalt_bak;rho = rho_bak; 
 deallocate(u_bak,v_bak,w_bak,du_bak,dv_bak,psi_bak,dpsi_bak)
 deallocate(temp_bak,salt_bak,dtemp_bak,dsalt_bak,rho_bak)
 
 first = .false.
 end subroutine diag_opt_balance
 
 
 
 
subroutine diag_opt_time_step
!=======================================================================
! rudimentary single time step
!=======================================================================
 use main_module   
 use diag_opt_balance_module
 implicit none
 integer :: i,j,k
 real*8 :: fxa, get_rho
 real*8,dimension(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx) :: fpx,fpy,forc

 do j=js_pe,je_pe
  do i=is_pe,ie_pe
   du_cor(i,j,:)= maskU(i,j,:)*(  coriolis_t(i  ,j)*(v(i  ,j,:,tau)+v(i  ,j-1,:,tau))*dxt(i  )/dxu(i)  &
                                 +coriolis_t(i+1,j)*(v(i+1,j,:,tau)+v(i+1,j-1,:,tau))*dxt(i+1)/dxu(i) )*0.25
   dv_cor(i,j,:)=-maskV(i,j,:)*(coriolis_t(i,j  )*(u(i-1,j  ,:,tau)+u(i,j  ,:,tau))*dyt(j  )*cost(j  )/( dyu(j)*cosu(j) )  &
                               +coriolis_t(i,j+1)*(u(i-1,j+1,:,tau)+u(i,j+1,:,tau))*dyt(j+1)*cost(j+1)/( dyu(j)*cosu(j) ) )*0.25
  enddo
 enddo

 call momentum_advection
 du(:,:,:,tau) = du_cor + opt_rho*du_adv
 dv(:,:,:,tau) = dv_cor + opt_rho*dv_adv

 !hydrostatic pressure
 fxa = grav/rho_0
 p_hydro(:,:,nz) = 0.5*rho(:,:,nz,tau)*fxa*dzw(nz)*maskT(:,:,nz)
 do k=nz-1,1,-1
   p_hydro(:,:,k)= maskT(:,:,k)*(p_hydro(:,:,k+1)+ 0.5*(rho(:,:,k+1,tau)+rho(:,:,k,tau))*fxa*dzw(k))
 enddo

 ! add hydrostatic pressure gradient to tendencies
 do j=js_pe,je_pe
  do i=is_pe,ie_pe
   du(i,j,:,tau) = du(i,j,:,tau) - ( p_hydro(i+1,j,:)-p_hydro(i,j,:)  )/(dxu(i)*cost(j)) *maskU(i,j,:) 
   dv(i,j,:,tau) = dv(i,j,:,tau) - ( p_hydro(i,j+1,:)-p_hydro(i,j,:)  ) /dyu(j)*maskV(i,j,:) 
  enddo
 enddo



 ! integrate forward in time
 u(:,:,:,taup1)=u(:,:,:,tau)+dt_mom*(  opt_A_ab*du(:,:,:,tau) + opt_B_ab*du(:,:,:,taum1) )*maskU
 v(:,:,:,taup1)=v(:,:,:,tau)+dt_mom*(  opt_A_ab*dv(:,:,:,tau) + opt_B_ab*dv(:,:,:,taum1) )*maskV
 
 if (enable_biharmonic_friction) then
   du_mix=0;dv_mix=0
   call biharmonic_friction
   u(:,:,:,taup1)=u(:,:,:,taup1)+dt_mom*du_mix
   v(:,:,:,taup1)=v(:,:,:,taup1)+dt_mom*dv_mix 
 endif 
 
 
 ! forcing for surface pressure 
 fpx=0.;fpy=0.
 do k=1,nz
    do j=js_pe,je_pe
      do i=is_pe,ie_pe
       fpx(i,j)=fpx(i,j)+u(i,j,k,taup1)*maskU(i,j,k)*dzt(k)/dt_mom
       fpy(i,j)=fpy(i,j)+v(i,j,k,taup1)*maskV(i,j,k)*dzt(k)/dt_mom
      enddo
     enddo
 enddo
 call border_exchg_xy(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,fpx) 
 call setcyclic_xy   (is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,fpx)
 call border_exchg_xy(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,fpy)
 call setcyclic_xy   (is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,fpy)

 ! forc = 1/cos (u_x + (cos v)_y )A
 do j=js_pe,je_pe
   do i=is_pe,ie_pe
    forc(i,j)=(fpx(i,j)-fpx(i-1,j))/(cost(j)*dxt(i))+(cosu(j)*fpy(i,j)-cosu(j-1)*fpy(i,j-1))/(cost(j)*dyt(j))
   enddo
 enddo
 if (enable_free_surface) then
  do j=js_pe,je_pe
   do i=is_pe,ie_pe
     forc(i,j) = forc(i,j) - psi(i,j,tau)/(grav*dt_mom**2)*maskT(i,j,nz)
   enddo
  enddo
 endif

 psi(:,:,taup1)=2*psi(:,:,tau)-psi(:,:,taum1) ! first guess
 !solve for surface pressure
 call congrad_surf_press(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,forc,congr_itts)
 call border_exchg_xy(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,psi(:,:,taup1)); 
 call setcyclic_xy(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,psi(:,:,taup1))

 ! remove surface pressure gradient
 do j=js_pe,je_pe
  do i=is_pe,ie_pe
   u(i,j,:,taup1) = u(i,j,:,taup1) - dt_mom*( psi(i+1,j,taup1)-psi(i,j,taup1))/(dxu(i)*cost(j)) *maskU(i,j,:) 
   v(i,j,:,taup1) = v(i,j,:,taup1) - dt_mom*( psi(i,j+1,taup1)-psi(i,j,taup1)) /dyu(j)*maskV(i,j,:) 
  enddo
 enddo
 call border_exchg_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,u(:,:,:,taup1)) 
 call setcyclic_xyz   (is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,u(:,:,:,taup1))
 call border_exchg_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,v(:,:,:,taup1)) 
 call setcyclic_xyz   (is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,v(:,:,:,taup1))
  
 ! integrate from bottom to surface to see error in w
 k=1
 do j=js_pe-onx+1,je_pe+onx
   do i=is_pe-onx+1,ie_pe+onx
         w(i,j,k,taup1) =-maskW(i,j,k)*dzt(k)* &
               ((        u(i,j,k,taup1)-          u(i-1,j,k,taup1))/(cost(j)*dxt(i)) &
               +(cosu(j)*v(i,j,k,taup1)-cosu(j-1)*v(i,j-1,k,taup1))/(cost(j)*dyt(j)) )
   enddo
 enddo
 do k=2,nz
  do j=js_pe-onx+1,je_pe+onx
   do i=is_pe-onx+1,ie_pe+onx
          w(i,j,k,taup1) = w(i,j,k-1,taup1)-maskW(i,j,k)*dzt(k)* &
               ((        u(i,j,k,taup1)          -u(i-1,j,k,taup1))/(cost(j)*dxt(i)) &
               +(cosu(j)*v(i,j,k,taup1)-cosu(j-1)*v(i,j-1,k,taup1))/(cost(j)*dyt(j)) )
   enddo
  enddo
 enddo

 
 call advect_tracer(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,temp_ave,dtemp_ave(:,:,:,tau) )
 call advect_tracer(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,temp(:,:,:,tau),dtemp(:,:,:,tau))
 temp(:,:,:,taup1)=temp(:,:,:,tau)+dt_tracer*( opt_A_ab*(dtemp_ave(:,:,:,tau)  +opt_rho*dtemp(:,:,:,tau)) &
                                             + opt_B_ab*(dtemp_ave(:,:,:,taum1)+opt_rho*dtemp(:,:,:,taum1)) )*maskT
 call border_exchg_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,temp(:,:,:,taup1)) 
 call setcyclic_xyz   (is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,temp(:,:,:,taup1))

 if (.not.opt_balance_temp_only) then
   call advect_tracer(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,salt_ave,dsalt_ave(:,:,:,tau))
   call advect_tracer(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,salt(:,:,:,tau),dsalt(:,:,:,tau))
   salt(:,:,:,taup1)=salt(:,:,:,tau)+dt_tracer*( opt_A_ab*(dsalt_ave(:,:,:,tau)  +opt_rho*dsalt(:,:,:,tau)) &
                                               + opt_B_ab*(dsalt_ave(:,:,:,taum1)+opt_rho*dsalt(:,:,:,taum1)))*maskT
   call border_exchg_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,salt(:,:,:,taup1)) 
   call setcyclic_xyz   (is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,salt(:,:,:,taup1))
   do k=1,nz
    do j=js_pe-onx,je_pe+onx
     do i=is_pe-onx,ie_pe+onx
      rho(i,j,k,taup1) = get_rho(salt(i,j,k,taup1)+salt_ave(i,j,k), &
                                 temp(i,j,k,taup1)+temp_ave(i,j,k),abs(zt(k)))*maskT(i,j,k)
    enddo
   enddo
  enddo
 else
   do k=1,nz
    do j=js_pe-onx,je_pe+onx
     do i=is_pe-onx,ie_pe+onx
      rho(i,j,k,taup1) = get_rho(35d0,temp(i,j,k,taup1)+temp_ave(i,j,k),abs(zt(k)))*maskT(i,j,k)
    enddo
   enddo
  enddo
 endif
 
end subroutine diag_opt_time_step





subroutine diag_opt_time_step_linear
!=======================================================================
! rudimentary single time step, but linear
!=======================================================================
 use main_module   
 use diag_opt_balance_module
 implicit none
 integer :: i,j,k
 real*8 :: fxa, get_rho
 real*8,dimension(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx) :: fpx,fpy,forc

 do j=js_pe,je_pe
  do i=is_pe,ie_pe
   du(i,j,:,tau)= maskU(i,j,:)*(  coriolis_t(i  ,j)*(v(i  ,j,:,tau)+v(i  ,j-1,:,tau))*dxt(i  )/dxu(i)  &
                                 +coriolis_t(i+1,j)*(v(i+1,j,:,tau)+v(i+1,j-1,:,tau))*dxt(i+1)/dxu(i) )*0.25
   dv(i,j,:,tau)=-maskV(i,j,:)*(coriolis_t(i,j  )*(u(i-1,j  ,:,tau)+u(i,j  ,:,tau))*dyt(j  )*cost(j  )/( dyu(j)*cosu(j) )  &
                               +coriolis_t(i,j+1)*(u(i-1,j+1,:,tau)+u(i,j+1,:,tau))*dyt(j+1)*cost(j+1)/( dyu(j)*cosu(j) ) )*0.25
  enddo
 enddo

 !hydrostatic pressure
 fxa = grav/rho_0
 p_hydro(:,:,nz) = 0.5*rho(:,:,nz,tau)*fxa*dzw(nz)*maskT(:,:,nz)
 do k=nz-1,1,-1
   p_hydro(:,:,k)= maskT(:,:,k)*(p_hydro(:,:,k+1)+ 0.5*(rho(:,:,k+1,tau)+rho(:,:,k,tau))*fxa*dzw(k))
 enddo

 ! add hydrostatic pressure gradient to tendencies
 do j=js_pe,je_pe
  do i=is_pe,ie_pe
   du(i,j,:,tau) = du(i,j,:,tau) - ( p_hydro(i+1,j,:)-p_hydro(i,j,:)  )/(dxu(i)*cost(j)) *maskU(i,j,:) 
   dv(i,j,:,tau) = dv(i,j,:,tau) - ( p_hydro(i,j+1,:)-p_hydro(i,j,:)  ) /dyu(j)*maskV(i,j,:) 
  enddo
 enddo

 ! integrate forward in time
 u(:,:,:,taup1)=u(:,:,:,tau)+dt_mom*(  opt_A_ab*du(:,:,:,tau) + opt_B_ab*du(:,:,:,taum1) )*maskU
 v(:,:,:,taup1)=v(:,:,:,tau)+dt_mom*(  opt_A_ab*dv(:,:,:,tau) + opt_B_ab*dv(:,:,:,taum1) )*maskV
 
 ! forcing for surface pressure 
 fpx=0.;fpy=0.
 do k=1,nz
    do j=js_pe,je_pe
      do i=is_pe,ie_pe
       fpx(i,j)=fpx(i,j)+u(i,j,k,taup1)*maskU(i,j,k)*dzt(k)/dt_mom
       fpy(i,j)=fpy(i,j)+v(i,j,k,taup1)*maskV(i,j,k)*dzt(k)/dt_mom
      enddo
     enddo
 enddo
 call border_exchg_xy(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,fpx) 
 call setcyclic_xy   (is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,fpx)
 call border_exchg_xy(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,fpy)
 call setcyclic_xy   (is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,fpy)

 ! forc = 1/cos (u_x + (cos v)_y )A
 do j=js_pe,je_pe
   do i=is_pe,ie_pe
    forc(i,j)=(fpx(i,j)-fpx(i-1,j))/(cost(j)*dxt(i))+(cosu(j)*fpy(i,j)-cosu(j-1)*fpy(i,j-1))/(cost(j)*dyt(j))
   enddo
 enddo
 if (enable_free_surface) then
  do j=js_pe,je_pe
   do i=is_pe,ie_pe
     forc(i,j) = forc(i,j) - psi(i,j,tau)/(grav*dt_mom**2)*maskT(i,j,nz)
   enddo
  enddo
 endif

 psi(:,:,taup1) = 2*psi(:,:,tau) - psi(:,:,taum1) ! first guess
 !solve for surface pressure
 call congrad_surf_press(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,forc,congr_itts)
 call border_exchg_xy(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,psi(:,:,taup1)); 
 call setcyclic_xy(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,psi(:,:,taup1))

 ! remove surface pressure gradient
 do j=js_pe,je_pe
  do i=is_pe,ie_pe
   u(i,j,:,taup1) = u(i,j,:,taup1) - dt_mom*( psi(i+1,j,taup1)-psi(i,j,taup1))/(dxu(i)*cost(j)) *maskU(i,j,:) 
   v(i,j,:,taup1) = v(i,j,:,taup1) - dt_mom*( psi(i,j+1,taup1)-psi(i,j,taup1)) /dyu(j)*maskV(i,j,:) 
  enddo
 enddo
 call border_exchg_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,u(:,:,:,taup1)) 
 call setcyclic_xyz   (is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,u(:,:,:,taup1))
 call border_exchg_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,v(:,:,:,taup1)) 
 call setcyclic_xyz   (is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,v(:,:,:,taup1))
  
 ! integrate from bottom to surface to see error in w
 k=1
 do j=js_pe-onx+1,je_pe+onx
   do i=is_pe-onx+1,ie_pe+onx
         w(i,j,k,taup1) =-maskW(i,j,k)*dzt(k)* &
               ((        u(i,j,k,taup1)-          u(i-1,j,k,taup1))/(cost(j)*dxt(i)) &
               +(cosu(j)*v(i,j,k,taup1)-cosu(j-1)*v(i,j-1,k,taup1))/(cost(j)*dyt(j)) )
   enddo
 enddo
 do k=2,nz
  do j=js_pe-onx+1,je_pe+onx
   do i=is_pe-onx+1,ie_pe+onx
          w(i,j,k,taup1) = w(i,j,k-1,taup1)-maskW(i,j,k)*dzt(k)* &
               ((        u(i,j,k,taup1)          -u(i-1,j,k,taup1))/(cost(j)*dxt(i)) &
               +(cosu(j)*v(i,j,k,taup1)-cosu(j-1)*v(i,j-1,k,taup1))/(cost(j)*dyt(j)) )
   enddo
  enddo
 enddo
 
 call advect_tracer(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,temp_ave,dtemp_ave(:,:,:,tau) )
 temp(:,:,:,taup1)=temp(:,:,:,tau)+dt_tracer*( opt_A_ab*(dtemp_ave(:,:,:,tau) ) &
                                             + opt_B_ab*(dtemp_ave(:,:,:,taum1)) )*maskT
 call border_exchg_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,temp(:,:,:,taup1)) 
 call setcyclic_xyz   (is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,temp(:,:,:,taup1))


 if (.not.opt_balance_temp_only) then
   call advect_tracer(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,salt_ave,dsalt_ave(:,:,:,tau))
   salt(:,:,:,taup1)=salt(:,:,:,tau)+dt_tracer*( opt_A_ab*(dsalt_ave(:,:,:,tau)  ) &
                                               + opt_B_ab*(dsalt_ave(:,:,:,taum1)))*maskT
   call border_exchg_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,salt(:,:,:,taup1)) 
   call setcyclic_xyz   (is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,salt(:,:,:,taup1))
   do k=1,nz
    do j=js_pe-onx,je_pe+onx
     do i=is_pe-onx,ie_pe+onx
      rho(i,j,k,taup1) = get_rho(salt(i,j,k,taup1)+salt_ave(i,j,k), &
                                 temp(i,j,k,taup1)+temp_ave(i,j,k),abs(zt(k)))*maskT(i,j,k)
     enddo
    enddo
   enddo
 else
   do k=1,nz
    do j=js_pe-onx,je_pe+onx
     do i=is_pe-onx,ie_pe+onx
      rho(i,j,k,taup1) = get_rho(35d0,temp(i,j,k,taup1)+temp_ave(i,j,k),abs(zt(k)))*maskT(i,j,k)
     enddo
    enddo
   enddo
 endif
end subroutine diag_opt_time_step_linear





subroutine diag_opt_time_ave(n_end,u_loc,v_loc,w_loc,temp_loc,salt_loc,rho_loc,psi_loc)
  use main_module
  use diag_opt_balance_module
  implicit none
  integer :: n,n_end
  real*8,dimension(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) :: u_loc,v_loc,w_loc,temp_loc,salt_loc,rho_loc
  real*8,dimension(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx)    :: psi_loc
 
  opt_A_ab = 1.; opt_B_ab = 0.;  taum1 = 1 ; tau = 2; taup1 = 3
  dt_tracer = + abs(dt_tracer); dt_mom = + abs(dt_mom)
  u_loc=0;v_loc=0;w_loc=0;temp_loc=0;salt_loc=0;rho_loc=0.;psi_loc=0.

  do n=0,n_end-1
    call diag_opt_time_step_linear
    opt_A_ab = 1.5 + AB_eps
    opt_B_ab = -(0.5 + AB_eps)   
    taup1   = mod(taup1,3)+1
    tau   = mod(tau,3)+1
    taum1 = mod(taum1,3)+1   
    u_loc = u_loc + u(:,:,:,tau)
    v_loc = v_loc + v(:,:,:,tau)
    w_loc = w_loc + w(:,:,:,tau)
    temp_loc = temp_loc + temp(:,:,:,tau)
    salt_loc = salt_loc + salt(:,:,:,tau)
    rho_loc  = rho_loc  + rho(:,:,:,tau)
    psi_loc  = psi_loc  + psi(:,:,tau)
  enddo
  
  u_loc = u_loc/n_end
  v_loc = v_loc/n_end
  w_loc = w_loc/n_end
  temp_loc = temp_loc/n_end
  salt_loc = salt_loc/n_end
  rho_loc  = rho_loc/n_end
  psi_loc  = psi_loc/n_end
  
  ! leave with well defined model state
  u(:,:,:,taup1) = u(:,:,:,tau); v(:,:,:,taup1) = v(:,:,:,tau); w(:,:,:,taup1) = w(:,:,:,tau)
  temp(:,:,:,taup1) = temp(:,:,:,tau); salt(:,:,:,taup1) = salt(:,:,:,tau);
  psi(:,:,taup1) = psi(:,:,tau); rho(:,:,:,taup1) = rho(:,:,:,tau)
  u(:,:,:,taum1) = u(:,:,:,tau); v(:,:,:,taum1) = v(:,:,:,tau); w(:,:,:,taum1) = w(:,:,:,tau)
  temp(:,:,:,taum1) = temp(:,:,:,tau); salt(:,:,:,taum1) = salt(:,:,:,tau);
  psi(:,:,taum1) = psi(:,:,tau); rho(:,:,:,taum1) = rho(:,:,:,tau)
  
end subroutine diag_opt_time_ave




subroutine diag_opt_backward_integration(n_end)
  use main_module
  use diag_opt_balance_module
  implicit none
  real*8 :: ramp,a_hbi_back
  integer :: n,n_end
  logical :: enable_back 
  
  enable_back = enable_superbee_advection 
  if (enable_superbee_advection ) then
   if (my_pe==0) print*,' switching off superbee advection for backward integration'
   enable_superbee_advection     = .false.
  endif
  if (a_hbi>0d0) then
   if (my_pe==0) print*,' switching off biharm. friction for backward integration'
   a_hbi_back = a_hbi
   a_hbi = 0d0
  endif
  
  opt_A_ab = 1.; opt_B_ab = 0.;  taum1 = 1 ; tau = 2; taup1 = 3
  dt_tracer = - abs(dt_tracer); dt_mom = - abs(dt_mom)
  do n=0,n_end-1
    opt_rho = 1. - ramp(n * dt_tracer /opt_balance_period)
    call diag_opt_time_step
    opt_A_ab = 1.5 + AB_eps
    opt_B_ab = -(0.5 + AB_eps)
    taup1   = mod(taup1,3)+1
    tau   = mod(tau,3)+1
    taum1 = mod(taum1,3)+1
  enddo
  
  ! leave with well defined model state
  u(:,:,:,taup1) = u(:,:,:,tau); v(:,:,:,taup1) = v(:,:,:,tau); w(:,:,:,taup1) = w(:,:,:,tau)
  temp(:,:,:,taup1) = temp(:,:,:,tau); salt(:,:,:,taup1) = salt(:,:,:,tau);
  psi(:,:,taup1) = psi(:,:,tau); rho(:,:,:,taup1) = rho(:,:,:,tau)
  u(:,:,:,taum1) = u(:,:,:,tau); v(:,:,:,taum1) = v(:,:,:,tau); w(:,:,:,taum1) = w(:,:,:,tau)
  temp(:,:,:,taum1) = temp(:,:,:,tau); salt(:,:,:,taum1) = salt(:,:,:,tau);
  psi(:,:,taum1) = psi(:,:,tau); rho(:,:,:,taum1) = rho(:,:,:,tau)
  
  enable_superbee_advection  = enable_back 
  a_hbi = a_hbi_back
end subroutine diag_opt_backward_integration




subroutine diag_opt_forward_integration(n_end)
  use main_module
  use diag_opt_balance_module
  implicit none
  real*8 :: ramp
  integer :: n,n_end
  opt_A_ab = 1.; opt_B_ab = 0.;  taum1 = 1 ; tau = 2; taup1 = 3
  dt_tracer = +abs(dt_tracer); dt_mom = + abs(dt_mom); 
  do n=0,n_end-1
    opt_rho = ramp(n * dt_tracer / opt_balance_period)
    call diag_opt_time_step
    opt_A_ab = 1.5 + AB_eps
    opt_B_ab = -(0.5 + AB_eps)
    taup1   = mod(taup1,3)+1
    tau   = mod(tau,3)+1
    taum1 = mod(taum1,3)+1
  enddo
  
  ! leave with well defined model state
  u(:,:,:,taup1) = u(:,:,:,tau); v(:,:,:,taup1) = v(:,:,:,tau); w(:,:,:,taup1) = w(:,:,:,tau)
  temp(:,:,:,taup1) = temp(:,:,:,tau); salt(:,:,:,taup1) = salt(:,:,:,tau);
  psi(:,:,taup1) = psi(:,:,tau); rho(:,:,:,taup1) = rho(:,:,:,tau)
  u(:,:,:,taum1) = u(:,:,:,tau); v(:,:,:,taum1) = v(:,:,:,tau); w(:,:,:,taum1) = w(:,:,:,tau)
  temp(:,:,:,taum1) = temp(:,:,:,tau); salt(:,:,:,taum1) = salt(:,:,:,tau);
  psi(:,:,taum1) = psi(:,:,tau); rho(:,:,:,taum1) = rho(:,:,:,tau)
  
end subroutine diag_opt_forward_integration



subroutine diag_opt_balance_norm(u_loc,v_loc)
  use main_module
  use diag_opt_balance_module
  implicit none
  real*8,dimension(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) :: u_loc,v_loc
  real*8 :: e1, e2, d
  integer :: i, j, k

  e1=0d0; e2=0d0; d=0d0
  do k=1,nz
    do j=js_pe,je_pe
      do i=is_pe,ie_pe
        e1 = e1 + u(i,j,k,tau)**2 + v(i,j,k,tau)**2 
        e2 = e2 + u_loc(i,j,k)**2 + v_loc(i,j,k)**2 
        d = d + (u(i,j,k,tau) - u_loc(i,j,k))**2 + (v(i,j,k,tau) - v_loc(i,j,k))**2           
      end do
    end do
  end do
  call global_sum(e1)
  call global_sum(e2)
  call global_sum(d) 
  norm_diff = 2d0 * sqrt(d) / (sqrt(e1) + sqrt(e2)) 
end subroutine diag_opt_balance_norm



subroutine diag_opt_write
 use main_module
 use diag_opt_balance_module
 implicit none
 include "netcdf.inc"
 integer :: ncid,iret,ilen,k,itimeid
 integer :: lon_tdim,lon_udim,itimedim,lat_tdim,lat_udim,id,z_tdim,z_udim
 real*8:: spval = -1.0d33,bloc(nx,ny),fxa
 logical,save :: first = .true.

 if (my_pe==0) print*,' writing to file opt_balance.cdf'

 if (first) then
  if (my_pe==0) print*,' preparing file opt_balance.cdf'
  call def_grid_cdf('opt_balance.cdf') 
  if (my_pe==0) then
    iret=nf_open('opt_balance.cdf',NF_WRITE, ncid)
    iret=nf_set_fill(ncid, NF_NOFILL, iret)
    call ncredf(ncid, iret)
    iret=nf_inq_dimid(ncid,'xt',lon_tdim)
    iret=nf_inq_dimid(ncid,'xu',lon_udim)
    iret=nf_inq_dimid(ncid,'yt',lat_tdim)
    iret=nf_inq_dimid(ncid,'yu',lat_udim)
    iret=nf_inq_dimid(ncid,'zt',z_tdim)
    iret=nf_inq_dimid(ncid,'zu',z_udim)
    iret=nf_inq_dimid(ncid,'Time',itimedim)
    
    id  = ncvdef (ncid,'psi',NCFLOAT,3,(/Lon_tdim,lat_tdim,iTimedim/),iret)
    call dvcdf(ncid,id,'surface pressure',len_trim('surface pressure'),'m^2/s^2',len_trim('m^2/s^2'),spval)
    id  = ncvdef (ncid,'u',NCFLOAT,4,(/Lon_udim,lat_tdim,z_tdim,iTimedim/),iret)
    call dvcdf(ncid,id,'zonal velocity',len_trim('zonal velocity'),'m/s',len_trim('m/s'),spval)
    id  = ncvdef (ncid,'v',NCFLOAT,4,(/Lon_tdim,lat_udim,z_tdim,iTimedim/),iret)
    call dvcdf(ncid,id,'meridional velocity',len_trim('meridional velocity'),'m/s',len_trim('m/s'),spval)
    id  = ncvdef (ncid,'w',NCFLOAT,4,(/Lon_tdim,lat_tdim,z_udim,iTimedim/),iret)
    call dvcdf(ncid,id,'vertical velocity',len_trim('vertical velocity'),'m/s',len_trim('m/s'),spval)
    id  = ncvdef (ncid,'temp',NCFLOAT,4,(/Lon_tdim,lat_tdim,z_tdim,iTimedim/),iret)
    call dvcdf(ncid,id,'temperature',len_trim('temperature'),'deg C',len_trim('deg C'),spval)
    id  = ncvdef (ncid,'salt',NCFLOAT,4,(/Lon_tdim,lat_tdim,z_tdim,iTimedim/),iret)
    call dvcdf(ncid,id,'salinity',len_trim('salinity'),'g/kg',len_trim('g/kg'),spval)

    id  = ncvdef (ncid,'psi_base',NCFLOAT,3,(/Lon_tdim,lat_tdim,iTimedim/),iret)
    call dvcdf(ncid,id,'surface pressure',len_trim('surface pressure'),'m^2/s^2',len_trim('m^2/s^2'),spval)
    id  = ncvdef (ncid,'u_base',NCFLOAT,4,(/Lon_udim,lat_tdim,z_tdim,iTimedim/),iret)
    call dvcdf(ncid,id,'zonal velocity',len_trim('zonal velocity'),'m/s',len_trim('m/s'),spval)
    id  = ncvdef (ncid,'v_base',NCFLOAT,4,(/Lon_tdim,lat_udim,z_tdim,iTimedim/),iret)
    call dvcdf(ncid,id,'meridional velocity',len_trim('meridional velocity'),'m/s',len_trim('m/s'),spval)
    id  = ncvdef (ncid,'w_base',NCFLOAT,4,(/Lon_tdim,lat_tdim,z_udim,iTimedim/),iret)
    call dvcdf(ncid,id,'vertical velocity',len_trim('vertical velocity'),'m/s',len_trim('m/s'),spval)
    id  = ncvdef (ncid,'temp_base',NCFLOAT,4,(/Lon_tdim,lat_tdim,z_tdim,iTimedim/),iret)
    call dvcdf(ncid,id,'temperature',len_trim('temperature'),'deg C',len_trim('deg C'),spval)
    id  = ncvdef (ncid,'salt_base',NCFLOAT,4,(/Lon_tdim,lat_tdim,z_tdim,iTimedim/),iret)
    call dvcdf(ncid,id,'salinity',len_trim('salinity'),'g/kg',len_trim('g/kg'),spval)
  
    id  = ncvdef (ncid,'temp_ave',NCFLOAT,1,(/z_tdim/),iret)
    call dvcdf(ncid,id,'temperature',len_trim('temperature'),'deg C',len_trim('deg C'),spval)
    id  = ncvdef (ncid,'salt_ave',NCFLOAT,1,(/z_tdim/),iret)
    call dvcdf(ncid,id,'salinity',len_trim('salinity'),'g/kg',len_trim('g/kg'),spval)
    
    call ncendf(ncid, iret)    
    iret=nf_inq_varid(ncid,'temp_ave',id)
    iret= nf_put_vara_double(ncid,id,(/1/), (/nz/),temp_ave(is_pe,js_pe,:))
    iret=nf_inq_varid(ncid,'salt_ave',id)
    iret= nf_put_vara_double(ncid,id,(/1/), (/nz/),salt_ave(is_pe,js_pe,:)) 
    iret=nf_close(ncid)     
  endif   
 endif   
 
 if (my_pe==0) then
    iret=nf_open('opt_balance.cdf',NF_WRITE, ncid)
    iret=nf_set_fill(ncid, NF_NOFILL, iret)
    iret=nf_inq_dimid(ncid,'Time',itimedim)
    iret=nf_inq_dimlen(ncid, itimedim,ilen)
    ilen = ilen+1
    fxa = itt*dt_tracer/86400.0
    if (fxa <1.0) then
     print'(a,f12.2,a,i4)',' writing snapshot at ',fxa*86400,' s, time steps in file : ',ilen
    else
     print'(a,f12.2,a,i4)',' writing snapshot at ',fxa,' d, time steps in file : ',ilen
    endif
    iret=nf_inq_varid(ncid,'Time',itimeid)
    iret= nf_put_vara_double(ncid,itimeid,(/ilen/),(/1/),(/fxa/))   
 endif    

 bloc(is_pe:ie_pe,js_pe:je_pe) = psi_bal(is_pe:ie_pe,js_pe:je_pe)
 where( maskT(is_pe:ie_pe,js_pe:je_pe,nz) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
 call pe0_recv_2D(nx,ny,bloc)
 if (my_pe==0) then
   iret=nf_inq_varid(ncid,'psi',id)
   iret= nf_put_vara_double(ncid,id,(/1,1,ilen/), (/nx,ny,1/),bloc)
 endif

 bloc(is_pe:ie_pe,js_pe:je_pe) = psi_base(is_pe:ie_pe,js_pe:je_pe)
 where( maskT(is_pe:ie_pe,js_pe:je_pe,nz) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
 call pe0_recv_2D(nx,ny,bloc)
 if (my_pe==0) then
   iret=nf_inq_varid(ncid,'psi_base',id)
   iret= nf_put_vara_double(ncid,id,(/1,1,ilen/), (/nx,ny,1/),bloc)
 endif
 
 do k=1,nz
     bloc(is_pe:ie_pe,js_pe:je_pe) = u_bal(is_pe:ie_pe,js_pe:je_pe,k)
     where( maskU(is_pe:ie_pe,js_pe:je_pe,k) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
     call pe0_recv_2D(nx,ny,bloc)
     if (my_pe==0) then
      iret=nf_inq_varid(ncid,'u',id)
      iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),bloc)
     endif
     bloc(is_pe:ie_pe,js_pe:je_pe) = u_base(is_pe:ie_pe,js_pe:je_pe,k)
     where( maskU(is_pe:ie_pe,js_pe:je_pe,k) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
     call pe0_recv_2D(nx,ny,bloc)
     if (my_pe==0) then
      iret=nf_inq_varid(ncid,'u_base',id)
      iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),bloc)
     endif
     bloc(is_pe:ie_pe,js_pe:je_pe) = v_bal(is_pe:ie_pe,js_pe:je_pe,k)
     where( maskV(is_pe:ie_pe,js_pe:je_pe,k) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
     call pe0_recv_2D(nx,ny,bloc)
     if (my_pe==0) then
      iret=nf_inq_varid(ncid,'v',id)
      iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),bloc)
     endif
     bloc(is_pe:ie_pe,js_pe:je_pe) = v_base(is_pe:ie_pe,js_pe:je_pe,k)
     where( maskV(is_pe:ie_pe,js_pe:je_pe,k) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
     call pe0_recv_2D(nx,ny,bloc)
     if (my_pe==0) then
      iret=nf_inq_varid(ncid,'v_base',id)
      iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),bloc)
     endif
     bloc(is_pe:ie_pe,js_pe:je_pe) = w_bal(is_pe:ie_pe,js_pe:je_pe,k)
     where( maskW(is_pe:ie_pe,js_pe:je_pe,k) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
     call pe0_recv_2D(nx,ny,bloc)
     if (my_pe==0) then
      iret=nf_inq_varid(ncid,'w',id)
      iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),bloc)
     endif
     bloc(is_pe:ie_pe,js_pe:je_pe) = w_base(is_pe:ie_pe,js_pe:je_pe,k)
     where( maskW(is_pe:ie_pe,js_pe:je_pe,k) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
     call pe0_recv_2D(nx,ny,bloc)
     if (my_pe==0) then
      iret=nf_inq_varid(ncid,'w_base',id)
      iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),bloc)
     endif
     bloc(is_pe:ie_pe,js_pe:je_pe) = temp_bal(is_pe:ie_pe,js_pe:je_pe,k)
     where( maskT(is_pe:ie_pe,js_pe:je_pe,k) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
     call pe0_recv_2D(nx,ny,bloc)
     if (my_pe==0) then
      iret=nf_inq_varid(ncid,'temp',id)
      iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),bloc)
     endif
     bloc(is_pe:ie_pe,js_pe:je_pe) = temp_base(is_pe:ie_pe,js_pe:je_pe,k)
     where( maskT(is_pe:ie_pe,js_pe:je_pe,k) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
     call pe0_recv_2D(nx,ny,bloc)
     if (my_pe==0) then
      iret=nf_inq_varid(ncid,'temp_base',id)
      iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),bloc)
     endif
     bloc(is_pe:ie_pe,js_pe:je_pe) = salt_bal(is_pe:ie_pe,js_pe:je_pe,k)
     where( maskT(is_pe:ie_pe,js_pe:je_pe,k) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
     call pe0_recv_2D(nx,ny,bloc)
     if (my_pe==0) then
      iret=nf_inq_varid(ncid,'salt',id)
      iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),bloc)
     endif
     bloc(is_pe:ie_pe,js_pe:je_pe) = salt_base(is_pe:ie_pe,js_pe:je_pe,k)
     where( maskT(is_pe:ie_pe,js_pe:je_pe,k) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
     call pe0_recv_2D(nx,ny,bloc)
     if (my_pe==0) then
      iret=nf_inq_varid(ncid,'salt_base',id)
      iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),bloc)
     endif
 enddo  
 if (my_pe==0) then
  iret=nf_close(ncid)
  print*,' done writing to file opt_balance.cdf'
 endif 
 
 first = .false.
end subroutine diag_opt_write 


function ramp(theta)
  implicit none 
  real*8, intent(in) :: theta
  real*8 :: ramp
  !real*8 :: pi = 3.14159265358979323846264338327950588
  real*8 :: t1,t2
  !ramp = (1d0-cos(pi*theta))*0.5d0
  t1 = 1./max(1d-32,theta )
  t2 = 1./max(1d-32,1d0-theta )
  ramp = exp(-t1)/(exp(-t1)+exp(-t2) )  
end function ramp





