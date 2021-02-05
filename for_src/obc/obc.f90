

subroutine set_obc_topo
 !---------------------------------------------------------------------------------
 ! routine is called in calc_topo
 !---------------------------------------------------------------------------------
 use main_module   
 use obc_module   
 implicit none
 integer :: i
 if (my_blk_j == 1 .and. enable_obc_south) then
   do i=is_pe-onx,ie_pe+onx
      kbot(i,1-onx:0)=kbot(i,1)
   enddo
 endif
 if (my_blk_j == n_pes_j .and. enable_obc_north) then
   do i=is_pe-onx,ie_pe+onx
    kbot(i,ny+1:ny+onx)=kbot(i,ny)
   enddo
 endif
end subroutine set_obc_topo


subroutine set_obc_boundary_xyz(is_,ie_,js_,je_,nz_,a)
 !---------------------------------------------------------------------------------
 ! routine is called in main, thermodynamics and calc_initial_conditions
 !---------------------------------------------------------------------------------
 use main_module   
 use obc_module   
 implicit none
 integer, intent(in) :: is_,ie_,js_,je_,nz_
 real*8, intent(inout)  :: a(is_:ie_,js_:je_,nz_)
 integer :: i,k
 if (my_blk_j == 1 .and. enable_obc_south) then
   do k=1,nz
    do i=is_pe-onx,ie_pe+onx
      a(i,1-onx:0,k)=a(i,1,k)
    enddo
   enddo
 endif
 if (my_blk_j == n_pes_j .and. enable_obc_north) then
   do k=1,nz
    do i=is_pe-onx,ie_pe+onx
      a(i,ny+1:ny+onx,k)=a(i,ny,k)
    enddo
   enddo
 endif
end subroutine set_obc_boundary_xyz



subroutine set_obc_streamfct
 !---------------------------------------------------------------------------------
 ! routine is called in solve_streamfunction
 !---------------------------------------------------------------------------------
 use main_module   
 use obc_module   
 implicit none
 integer :: i,j
 real*8 :: c1ps,c1pn, var, var1
 
 if (my_blk_j == 1 .and. enable_obc_south .and. .not. enable_prescribe_psi_obc_south ) then
   j=1
   var=-dyu(j+1)/dt_mom
   do i=is_pe,ie_pe
    var1=psi(i,j+2,tau)-psi(i,j+1,tau)
    if (var1== 0.) then
      c1ps=var
    else
      c1ps=var*(psi(i,j+1,taup1)-psi(i,j+1,tau))/var1
      if (c1ps.ge. 0.) c1ps=0.
      if (c1ps.lt.var) c1ps=var
    endif
    psi(i,j  ,taup1) = psi(i,j,tau)  - c1ps*dt_mom/dyu(j)*(psi(i,j+1,tau)-psi(i,j,tau))*maskZ(i,j,nz)
    psi(i,j-1,taup1) = psi(i,j,taup1)
   enddo
 endif
      
 if (my_blk_j == n_pes_j .and. enable_obc_north .and. .not. enable_prescribe_psi_obc_north) then
   j=ny
   var=dyu(j-1)/dt_mom
   do i=is_pe,ie_pe
    var1=psi(i,j-1,tau)-psi(i,j-2,tau)
    if (var1.eq.0.) then
      c1pn=var
    else
      c1pn=-var*(psi(i,j-1,taup1)-psi(i,j-1,tau))/var1
      if (c1pn.le. 0.) c1pn=0.
      if (c1pn.gt.var) c1pn=var
    endif
    psi(i,j  ,taup1)  = psi(i,j,tau)- c1pn*dt_mom/dyu(j)*(psi(i,j,tau)-psi(i,j-1,tau))*maskZ(i,j,nz)
    psi(i,j+1,taup1)  = psi(i,j,taup1)
   enddo
 endif
 
 if (my_blk_j == 1 .and. enable_obc_south .and. enable_prescribe_psi_obc_south) then
   do i=is_pe-onx,ie_pe+onx
     psi(i,1-onx:1,taup1) = psi_wall_south(i)*maskZ(i,1,nz)
   enddo
 endif

 if (my_blk_j == n_pes_j .and. enable_obc_north .and. enable_prescribe_psi_obc_north) then
   do i=is_pe-onx,ie_pe+onx
     psi(i,ny:ny+onx,taup1) = psi_wall_north(i)*maskZ(i,ny,nz)
   enddo
 endif

 if (enable_obc_south .or. enable_obc_north) then
  call border_exchg_xy(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,psi(:,:,taup1))
  call setcyclic_xy   (is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,psi(:,:,taup1))
 endif
 
end subroutine set_obc_streamfct




subroutine set_obc_momentum
 !---------------------------------------------------------------------------------
 ! routine is called in momentum
 !---------------------------------------------------------------------------------
 use main_module   
 use obc_module   
 implicit none
 integer :: i,j,k
 real*8 :: diss(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz),obc_ah(js_pe-onx:je_pe+onx)
 real*8 :: yN,yS
 
 !---------------------------------------------------------------------------------
 ! Linear momentum equation at open boundary
 !---------------------------------------------------------------------------------
 if (my_blk_j == 1       .and. enable_obc_south) then
    j=1
    du(:,j,:,tau) = du(:,j,:,tau) - du_adv(:,j,:)
    dv(:,j,:,tau) = dv(:,j,:,tau) - dv_adv(:,j,:)
 endif
 if (my_blk_j == n_pes_j .and. enable_obc_north) then
    j=ny
    du(:,j,:,tau) = du(:,j,:,tau) - du_adv(:,j,:)
    dv(:,j,:,tau) = dv(:,j,:,tau) - dv_adv(:,j,:)
 endif

 obc_ah=0.0
 
 if (enable_obc_north_damping) then
    yN = 0.0; if (my_blk_j == n_pes_j) yN = yu(ny)
    call global_max(yN)
    do j=js_pe-onx,je_pe+onx
       obc_ah(j) = obc_ah(j) + obc_north_damping_amp*exp( -(yN-yt(j))**2/obc_north_damping_len**2  )
    enddo
 endif
 
 if (enable_obc_south_damping) then
    yS = 1e20; if (my_blk_j == 1) yS = yu(0)
    call global_min(yS)
    do j=js_pe-onx,je_pe+onx
       obc_ah(j) = obc_ah(j) + obc_south_damping_amp*exp( -(yt(j)-yS)**2/obc_south_damping_len**2  )
    enddo
 endif

 !---------------------------------------------------------------------------------
 !  Damping zones with enhanced harmonic friction
 !---------------------------------------------------------------------------------
 if (enable_obc_south_damping.or.enable_obc_north_damping) then

  !---------------------------------------------------------------------------------
  ! zonal velocity
  !---------------------------------------------------------------------------------
  do j=js_pe,je_pe
   do i=is_pe-2,ie_pe
    flux_east(i,j,:)=obc_ah(j)*(u(i+1,j,:,tau)-u(i,j,:,tau))/(cost(j)*dxt(i+1))*maskU(i+1,j,:)*maskU(i,j,:)
   enddo
  enddo
  do j=js_pe-1,je_pe
   flux_north(:,j,:)=(obc_ah(j+1)+obc_ah(j))/2.0* &
             (u(:,j+1,:,tau)-u(:,j,:,tau))/dyu(j)*maskU(:,j+1,:)*maskU(:,j,:)*cosu(j)
  enddo 
  do j=js_pe,je_pe
   do i=is_pe,ie_pe
    du_mix(i,j,:)= du_mix(i,j,:) + maskU(i,j,:)*((flux_east(i,j,:) - flux_east(i-1,j,:))/(cost(j)*dxu(i)) &
                                                +(flux_north(i,j,:) - flux_north(i,j-1,:))/(cost(j)*dyt(j)) )
   enddo
  enddo
  if (enable_conserve_energy) then
   do k=1,nz
    do j=js_pe,je_pe
     do i=is_pe-1,ie_pe
      diss(i,j,k) =0.5*((u(i+1,j,k,tau)-u(i,j,k,tau))*flux_east(i,j,k) &
                       +(u(i,j,k,tau)-u(i-1,j,k,tau))*flux_east(i-1,j,k))/(cost(j)*dxu(i))  &
                  +0.5*((u(i,j+1,k,tau)-u(i,j,k,tau))*flux_north(i,j,k)+ &
                        (u(i,j,k,tau)-u(i,j-1,k,tau))*flux_north(i,j-1,k))/(cost(j)*dyt(j)) 
     enddo
    enddo
   enddo
   call calc_diss(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,diss,K_diss_h,'U')
  endif

  !---------------------------------------------------------------------------------
  ! meridional velocity
  !---------------------------------------------------------------------------------
  do j=js_pe-1,je_pe
   do i=is_pe-1,ie_pe
      flux_east(i,j,:)=(obc_ah(j+1)+obc_ah(j))/2.0* &
                (v(i+1,j,:,tau)-v(i,j,:,tau))/(cosu(j)*dxu(i)) *maskV(i+1,j,:)*maskV(i,j,:)
   enddo
  enddo
  do j=js_pe-2,je_pe
   flux_north(:,j,:)=obc_ah(j+1)*(v(:,j+1,:,tau)-v(:,j,:,tau) )/dyt(j+1)*cost(j+1)*maskV(:,j,:)*maskV(:,j+1,:)
  enddo
  do j=js_pe,je_pe
   do i=is_pe,ie_pe
    dv_mix(i,j,:)= dv_mix(i,j,:) + maskV(i,j,:)*( (flux_east(i,j,:) - flux_east(i-1,j,:))/(cosu(j)*dxt(i))  &
                                                 +(flux_north(i,j,:) - flux_north(i,j-1,:))/(dyu(j)*cosu(j)) )
   enddo
  enddo
  if (enable_conserve_energy) then
   do k=1,nz
    do j=js_pe-1,je_pe
     do i=is_pe,ie_pe
      diss(i,j,k) =0.5*((v(i+1,j,k,tau)-v(i,j,k,tau))*flux_east(i,j,k)+ &
                        (v(i,j,k,tau)-v(i-1,j,k,tau))*flux_east(i-1,j,k))/(cosu(j)*dxt(i)) &
                 + 0.5*((v(i,j+1,k,tau)-v(i,j,k,tau))*flux_north(i,j,k)+ &
                        (v(i,j,k,tau)-v(i,j-1,k,tau))*flux_north(i,j-1,k))/(cosu(j)*dyu(j)) 
     enddo
    enddo
   enddo
   call calc_diss(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,diss,K_diss_h,'V')
  endif
 endif
end subroutine set_obc_momentum

  
subroutine set_obc_tempsalt
 !---------------------------------------------------------------------------------
 ! routine is called in thermodynamic
 !---------------------------------------------------------------------------------
 use main_module   
 use obc_module   
 implicit none
 integer :: j
 real*8 :: var
 real*8, dimension(is_pe-onx:ie_pe+onx,nz) :: v_adv,var1,ph_vel
 
!-----------------------------------------------------------------------
!    1) compute the advective velocity "vad" at the north/south face of the "t" grid box
!    2) compute phase velocity at the boundary
!    3) calculate radiation condition at the wall
!-----------------------------------------------------------------------

 if (my_blk_j == 1 .and. enable_obc_south) then
    
    j=1
    var = -dyu(j+1)/dt_tracer
    v_adv = min( 0D0, v(:,j,:,tau)*maskT(:,j,:) )
    
    var1=temp(:,j+2,:,taum1)-temp(:,j+1,:,taum1)
    where ( var1 == 0d0) 
      ph_vel=var
    elsewhere
       ph_vel = var*(temp(:,j+1,:,tau)-temp(:,j+1,:,taum1))/var1*maskT(:,j,:)
       ph_vel = max( var, min(0D0, ph_vel)  )
    end where
    temp(:,j,:,taup1) = temp(:,j,:,tau) + dt_tracer*maskT(:,j,:)*   &
         (  -(ph_vel+v_adv)*(temp(:,j+1,:,tau)-temp(:,j,:,tau))/dyu(j) )

    var1=salt(:,j+2,:,taum1)-salt(:,j+1,:,taum1)
    where ( var1 == 0d0) 
      ph_vel=var
    elsewhere
       ph_vel = var*(salt(:,j+1,:,tau)-salt(:,j+1,:,taum1))/var1*maskT(:,j,:)
       ph_vel = max( var, min(0D0, ph_vel)  )
    end where
    salt(:,j,:,taup1) = salt(:,j,:,tau) + dt_tracer*maskT(:,j,:)*   &
         (  -(ph_vel+v_adv)*(salt(:,j+1,:,tau)-salt(:,j,:,tau))/dyu(j) )
    
   if (enable_restore_TS_obc_south) then
     temp(:,j,:,taup1) = temp(:,j,:,taup1) + dt_tracer*maskT(:,j,:)*obc_tscl*(temp_wall_south - temp(:,j,:,tau) )
     salt(:,j,:,taup1) = salt(:,j,:,taup1) + dt_tracer*maskT(:,j,:)*obc_tscl*(salt_wall_south - salt(:,j,:,tau) )
   endif
   call obc_hor_mixing(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,temp,j) 
   call obc_hor_mixing(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,salt,j) 
   
 endif

 if (my_blk_j == n_pes_j .and. enable_obc_north) then
   j=ny
   var=dyu(j-1)/dt_tracer
   v_adv = max( 0D0, v(:,j-1,:,tau)*maskT(:,j,:) )
   
   var1 = temp(:,j-1,:,taum1)-temp(:,j-2,:,taum1)
   where (var1.eq.0.) 
      ph_vel = var
   elsewhere
      ph_vel = -var*(temp(:,j-1,:,tau)-temp(:,j-1,:,taum1)) /var1*maskT(:,j,:)
      ph_vel = min(var, max( 0D0, ph_vel)  )
   end where
   temp(:,j,:,taup1) = temp(:,j,:,tau) + dt_tracer*maskT(:,j,:)*   &
           ( -(ph_vel+v_adv)*(temp(:,j,:,tau)-temp(:,j-1,:,tau))/dyu(j-1) )

   var1 = salt(:,j-1,:,taum1)-salt(:,j-2,:,taum1)
   where (var1.eq.0.) 
      ph_vel = var
   elsewhere
      ph_vel = -var*(salt(:,j-1,:,tau)-salt(:,j-1,:,taum1)) /var1*maskT(:,j,:)
      ph_vel = min(var, max( 0D0, ph_vel)  )
   end where
   salt(:,j,:,taup1) = salt(:,j,:,tau) + dt_tracer*maskT(:,j,:)*   &
           ( -(ph_vel+v_adv)*(salt(:,j,:,tau)-salt(:,j-1,:,tau))/dyu(j-1) )

   if (enable_restore_TS_obc_north) then
     temp(:,j,:,taup1) = temp(:,j,:,taup1) + dt_tracer*maskT(:,j,:)*obc_tscl*(temp_wall_north - temp(:,j,:,tau) )
     salt(:,j,:,taup1) = salt(:,j,:,taup1) + dt_tracer*maskT(:,j,:)*obc_tscl*(salt_wall_north - salt(:,j,:,tau) )
   endif
   call obc_hor_mixing(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,temp,j) 
   call obc_hor_mixing(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,salt,j) 
 endif       
end subroutine set_obc_tempsalt


subroutine set_obc_tracer(is_,ie_,js_,je_,nz_,a,positive_definite)
 !---------------------------------------------------------------------------------
 ! routine is called in ...
 !---------------------------------------------------------------------------------
 use main_module   
 use obc_module   
 implicit none
 integer, intent(in) :: is_,ie_,js_,je_,nz_
 real*8, intent(inout)  :: a(is_:ie_,js_:je_,nz_,3)
 integer :: j
 logical :: positive_definite
 real*8 :: var
 real*8, dimension(is_pe-onx:ie_pe+onx,nz) :: v_adv,var1,ph_vel
 
 if (my_blk_j == 1 .and. enable_obc_south) then
    j=1
    var = -dyu(j+1)/dt_tracer
    v_adv = min( 0D0, v(:,j,:,tau)*maskT(:,j,:) )
    var1=a(:,j+2,:,taum1)-a(:,j+1,:,taum1)
    where ( var1 == 0d0) 
      ph_vel=var
    elsewhere
       ph_vel = var*(a(:,j+1,:,tau)-a(:,j+1,:,taum1))/var1*maskT(:,j,:)
       ph_vel = max( var, min(0D0, ph_vel)  )
    end where
    a(:,j,:,taup1) = a(:,j,:,tau) + dt_tracer*maskT(:,j,:)*   &
         (  -(ph_vel+v_adv)*(a(:,j+1,:,tau)-a(:,j,:,tau))/dyu(j) )
    if (positive_definite) a(:,j,:,taup1) = max( 0d0,  a(:,j,:,taup1))
    call obc_hor_mixing(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,a,j) 
 endif

 if (my_blk_j == n_pes_j .and. enable_obc_north) then
   j=ny
   var=dyu(j-1)/dt_tracer
   v_adv = max( 0D0, v(:,j-1,:,tau)*maskT(:,j,:) )
   var1 = a(:,j-1,:,taum1)-a(:,j-2,:,taum1)
   where (var1.eq.0.) 
      ph_vel = var
   elsewhere
      ph_vel = -var*(a(:,j-1,:,tau)-a(:,j-1,:,taum1)) /var1*maskT(:,j,:)
      ph_vel = min(var, max( 0D0, ph_vel)  )
   end where
   a(:,j,:,taup1) = a(:,j,:,tau) + dt_tracer*maskT(:,j,:)*   &
           ( -(ph_vel+v_adv)*(a(:,j,:,tau)-a(:,j-1,:,tau))/dyu(j-1) )
   if (positive_definite) a(:,j,:,taup1) = max( 0d0,  a(:,j,:,taup1))
   call obc_hor_mixing(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,a,j) 
 endif       
end subroutine set_obc_tracer




subroutine set_obc_tracer_wgrid(is_,ie_,js_,je_,nz_,a,positive_definite)
 !---------------------------------------------------------------------------------
 ! routine is called in ...
 !---------------------------------------------------------------------------------
 use main_module   
 use obc_module   
 implicit none
 integer, intent(in) :: is_,ie_,js_,je_,nz_
 real*8, intent(inout)  :: a(is_:ie_,js_:je_,nz_,3)
 integer :: j
 logical :: positive_definite
 real*8 :: var
 real*8, dimension(is_pe-onx:ie_pe+onx,nz) :: v_adv,var1,ph_vel
 
 if (my_blk_j == 1 .and. enable_obc_south) then
    j=1
    var = -dyu(j+1)/dt_tracer
    v_adv = min( 0D0, v_wgrid(:,j,:)*maskW(:,j,:) )
    var1=a(:,j+2,:,taum1)-a(:,j+1,:,taum1)
    where ( var1 == 0d0) 
      ph_vel=var
    elsewhere
       ph_vel = var*(a(:,j+1,:,tau)-a(:,j+1,:,taum1))/var1*maskW(:,j,:)
       ph_vel = max( var, min(0D0, ph_vel)  )
    end where
    a(:,j,:,taup1) = a(:,j,:,tau) + dt_tracer*maskW(:,j,:)*   &
         (  -(ph_vel+v_adv)*(a(:,j+1,:,tau)-a(:,j,:,tau))/dyu(j) )
    if (positive_definite) a(:,j,:,taup1) = max( 0d0,  a(:,j,:,taup1))
    call obc_hor_mixing_wgrid(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,a,j) 
 endif

 if (my_blk_j == n_pes_j .and. enable_obc_north) then
   j=ny
   var=dyu(j-1)/dt_tracer
   v_adv = max( 0D0, v_wgrid(:,j-1,:)*maskW(:,j,:) )
   var1 = a(:,j-1,:,taum1)-a(:,j-2,:,taum1)
   where (var1.eq.0.) 
      ph_vel = var
   elsewhere
      ph_vel = -var*(a(:,j-1,:,tau)-a(:,j-1,:,taum1)) /var1*maskW(:,j,:)
      ph_vel = min(var, max( 0D0, ph_vel)  )
   end where
   a(:,j,:,taup1) = a(:,j,:,tau) + dt_tracer*maskW(:,j,:)*   &
           ( -(ph_vel+v_adv)*(a(:,j,:,tau)-a(:,j-1,:,tau))/dyu(j-1) )
   if (positive_definite) a(:,j,:,taup1) = max( 0d0,  a(:,j,:,taup1))
   call obc_hor_mixing_wgrid(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,a,j) 
 endif       
end subroutine set_obc_tracer_wgrid






subroutine obc_hor_mixing(is_,ie_,js_,je_,nz_,a,j)
 use main_module   
 use obc_module   
 implicit none
 integer, intent(in) :: is_,ie_,js_,je_,nz_,j
 real*8, intent(inout)  :: a(is_:ie_,js_:je_,nz_,3)
 integer :: i,j1
 do i=is_pe-1,ie_pe
     flux_east(i,j,:)=obc_K_h*(a(i+1,j,:,tau)-a(i,j,:,tau))/(cost(j)*dxu(i))*maskU(i,j,:)
 enddo
 do j1=j-1,j
   flux_north(:,j1,:)=obc_K_h*(a(:,j1+1,:,tau)-a(:,j1,:,tau))/dyu(j1)*maskV(:,j1,:)*cosu(j1)
 enddo
 do i=is_pe,ie_pe
    a(i,j,:,taup1)=a(i,j,:,taup1)+dt_tracer*maskT(i,j,:)  &
                                   *((flux_east(i,j,:) - flux_east(i-1,j,:))/(cost(j)*dxt(i)) &
                                    +(flux_north(i,j,:) -flux_north(i,j-1,:))/(cost(j)*dyt(j)) )
 enddo
end subroutine obc_hor_mixing





subroutine obc_hor_mixing_wgrid(is_,ie_,js_,je_,nz_,a,j)
 use main_module   
 use obc_module   
 implicit none
 integer, intent(in) :: is_,ie_,js_,je_,nz_,j
 real*8, intent(inout)  :: a(is_:ie_,js_:je_,nz_,3)
 integer :: i,j1,k
 real*8 :: maskUtr,maskVtr
 maskUtr(i,j,k) = maskW(i+1,j,k)*maskW(i,j,k)
 maskVtr(i,j,k) = maskW(i,j+1,k)*maskW(i,j,k)
 do k=1,nz
  do i=is_pe-1,ie_pe
     flux_east(i,j,k)=obc_K_h*(a(i+1,j,k,tau)-a(i,j,k,tau))/(cost(j)*dxu(i))*maskUtr(i,j,k)
  enddo
  do j1=j-1,j
   do i=is_pe,ie_pe
     flux_north(i,j1,k)=obc_K_h*(a(i,j1+1,k,tau)-a(i,j1,k,tau))/dyu(j1)*maskVtr(i,j1,k)*cosu(j1)
   enddo
  enddo
 enddo
 do i=is_pe,ie_pe
    a(i,j,:,taup1)=a(i,j,:,taup1)+dt_tracer*maskT(i,j,:)  &
                                   *((flux_east(i,j,:) - flux_east(i-1,j,:))/(cost(j)*dxt(i)) &
                                    +(flux_north(i,j,:) -flux_north(i,j-1,:))/(cost(j)*dyt(j)) )
 enddo
end subroutine obc_hor_mixing_wgrid
