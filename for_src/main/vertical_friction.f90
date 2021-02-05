



subroutine explicit_vert_friction
!=======================================================================
!  explicit vertical friction
!  dissipation is calculated and added to K_diss_v
!=======================================================================
 use main_module   
 implicit none
 integer :: i,j,k
 real*8 :: diss(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz),fxa
 
 !---------------------------------------------------------------------------------
 ! vertical friction of zonal momentum
 !---------------------------------------------------------------------------------
 do k=1,nz-1
  do j=js_pe-1,je_pe
   do i=is_pe-1,ie_pe
    fxa = 0.5*(kappaM(i,j,k)+kappaM(i+1,j,k))
    flux_top(i,j,k)=fxa*(u(i,j,k+1,tau)-u(i,j,k,tau))/dzw(k)*maskU(i,j,k+1)*maskU(i,j,k)
   enddo
  enddo
 enddo
 flux_top(:,:,nz)=0d0
 k=1; du_mix(:,:,k) = flux_top(:,:,k)/dzt(k)*maskU(:,:,k)
 do k=2,nz
   du_mix(:,:,k) = (flux_top(:,:,k)-flux_top(:,:,k-1))/dzt(k)*maskU(:,:,k)
 enddo

 !---------------------------------------------------------------------------------
 ! diagnose dissipation by vertical friction of zonal momentum 
 !---------------------------------------------------------------------------------
  do k=1,nz-1
   do j=js_pe-1,je_pe
    do i=is_pe-1,ie_pe
     diss(i,j,k) = (u(i,j,k+1,tau)-u(i,j,k,tau))*flux_top(i,j,k)/dzw(k)  
    enddo
   enddo
  enddo
  diss(:,:,nz)=0.0
  call ugrid_to_tgrid(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,diss,diss)
  K_diss_v = K_diss_v + diss 

 !---------------------------------------------------------------------------------
 ! vertical friction of meridional momentum
 !---------------------------------------------------------------------------------
 do k=1,nz-1
  do j=js_pe-1,je_pe
   do i=is_pe-1,ie_pe
    fxa = 0.5*(kappaM(i,j,k)+kappaM(i,j+1,k))
    flux_top(i,j,k)=fxa*(v(i,j,k+1,tau)-v(i,j,k,tau))/dzw(k)*maskV(i,j,k+1)*maskV(i,j,k)
   enddo
  enddo
 enddo
 flux_top(:,:,nz)=0d0
 k=1; dv_mix(:,:,k) = flux_top(:,:,k)/dzt(k)*maskV(:,:,k)
 do k=2,nz
   dv_mix(:,:,k) = (flux_top(:,:,k)-flux_top(:,:,k-1))/dzt(k)*maskV(:,:,k)
 enddo

 !---------------------------------------------------------------------------------
 ! diagnose dissipation by vertical friction of meridional momentum
 !---------------------------------------------------------------------------------
  do k=1,nz-1
   do j=js_pe-1,je_pe
    do i=is_pe-1,ie_pe
     diss(i,j,k) = (v(i,j,k+1,tau)-v(i,j,k,tau))*flux_top(i,j,k)/dzw(k)  
    enddo
   enddo
  enddo
  diss(:,:,nz)=0.0
  call vgrid_to_tgrid(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,diss,diss)
  K_diss_v = K_diss_v + diss 

 if (.not.enable_hydrostatic) then
 !---------------------------------------------------------------------------------
 ! vertical friction of vertical momentum
 !---------------------------------------------------------------------------------
   do k=1,nz-1
    do j=js_pe-1,je_pe
     do i=is_pe-1,ie_pe
      fxa = 0.5*(kappaM(i,j,k)+kappaM(i,j,k+1))
      flux_top(i,j,k)=fxa*(w(i,j,k+1,tau)-w(i,j,k,tau))/dzt(k+1)*maskW(i,j,k+1)*maskW(i,j,k)
     enddo
    enddo
   enddo
   flux_top(:,:,nz)=0d0
   k=1; dw_mix(:,:,k) = flux_top(:,:,k)/dzw(k)*maskW(:,:,k)
   do k=2,nz
    dw_mix(:,:,k) = (flux_top(:,:,k)-flux_top(:,:,k-1))/dzw(k)*maskW(:,:,k)
   enddo

 !---------------------------------------------------------------------------------
 ! diagnose dissipation by vertical friction of vertical momentum
 !---------------------------------------------------------------------------------
 ! to be implemented
 endif

end subroutine explicit_vert_friction






subroutine implicit_vert_friction
!=======================================================================
!  vertical friction
!  dissipation is calculated and added to K_diss_v
!=======================================================================
 use main_module   
 implicit none
 integer :: i,j,k,ks
 real*8 :: a_tri(nz),b_tri(nz),c_tri(nz),d_tri(nz),delta(nz),fxa
 real*8 :: diss(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz)
 
 !---------------------------------------------------------------------------------
 ! implicit vertical friction of zonal momentum
 !---------------------------------------------------------------------------------
 do j=js_pe-1,je_pe
   do i=is_pe-1,ie_pe
    ks=max(kbot(i,j),kbot(i+1,j))
    if (ks>0) then
     do k=ks,nz-1
      fxa = 0.5*(kappaM(i,j,k)+kappaM(i+1,j,k))
      delta(k) = dt_mom/dzw(k)*fxa*maskU(i,j,k+1)*maskU(i,j,k)
     enddo
     delta(nz)=0.0
     a_tri(ks)=0.0
     do k=ks+1,nz
       a_tri(k) = - delta(k-1)/dzt(k)
     enddo
     b_tri(ks) = 1+ delta(ks)/dzt(ks)   
     do k=ks+1,nz-1
      b_tri(k) = 1+ delta(k)/dzt(k) + delta(k-1)/dzt(k) 
     enddo
     b_tri(nz) = 1+ delta(nz-1)/dzt(nz) 
     do k=ks,nz-1
      c_tri(k) = - delta(k)/dzt(k)
     enddo
     c_tri(nz)=0.0
     d_tri(ks:nz)=u(i,j,ks:nz,tau) 
     call solve_tridiag(a_tri(ks:nz),b_tri(ks:nz),c_tri(ks:nz),d_tri(ks:nz),u(i,j,ks:nz,taup1),nz-ks+1)
    endif
    du_mix(i,j,:)=(u(i,j,:,taup1)-u(i,j,:,tau))/dt_mom
   enddo
 enddo

 !---------------------------------------------------------------------------------
 ! diagnose dissipation by vertical friction of zonal momentum 
 !---------------------------------------------------------------------------------
  do k=1,nz-1
   do j=js_pe-1,je_pe
    do i=is_pe-1,ie_pe
     fxa = 0.5*(kappaM(i,j,k)+kappaM(i+1,j,k))
     flux_top(i,j,k)=fxa*(u(i,j,k+1,taup1)-u(i,j,k,taup1))/dzw(k)*maskU(i,j,k+1)*maskU(i,j,k)
    enddo
   enddo
  enddo
  do k=1,nz-1
   do j=js_pe-1,je_pe
    do i=is_pe-1,ie_pe
     diss(i,j,k) = (u(i  ,j,k+1,tau)-u(i  ,j,k,tau))*flux_top(i  ,j,k)/dzw(k)  
    enddo
   enddo
  enddo
  diss(:,:,nz)=0.0
  call ugrid_to_tgrid(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,diss,diss)
  K_diss_v = K_diss_v + diss 

 !---------------------------------------------------------------------------------
 ! implicit vertical friction of meridional momentum
 !---------------------------------------------------------------------------------
 do j=js_pe-1,je_pe
   do i=is_pe-1,ie_pe
    ks=max(kbot(i,j),kbot(i,j+1))
    if (ks>0) then
     do k=ks,nz-1
      fxa = 0.5*(kappaM(i,j,k)+kappaM(i,j+1,k))
      delta(k) = dt_mom/dzw(k)*fxa*maskV(i,j,k+1)*maskV(i,j,k)
     enddo
     delta(nz)=0.0
     a_tri(ks)=0.0
     do k=ks+1,nz
       a_tri(k) = - delta(k-1)/dzt(k)
     enddo
     b_tri(ks) = 1+ delta(ks)/dzt(ks)   
     do k=ks+1,nz-1
      b_tri(k) = 1+ delta(k)/dzt(k) + delta(k-1)/dzt(k) 
     enddo
     b_tri(nz) = 1+ delta(nz-1)/dzt(nz) 
     do k=ks,nz-1
      c_tri(k) = - delta(k)/dzt(k)
     enddo
     c_tri(nz)=0.0
     d_tri(ks:nz)=v(i,j,ks:nz,tau) 
     call solve_tridiag(a_tri(ks:nz),b_tri(ks:nz),c_tri(ks:nz),d_tri(ks:nz),v(i,j,ks:nz,taup1),nz-ks+1)
    endif
    dv_mix(i,j,:)=(v(i,j,:,taup1)-v(i,j,:,tau))/dt_mom
   enddo
 enddo

 !---------------------------------------------------------------------------------
 ! diagnose dissipation by vertical friction of meridional momentum
 !---------------------------------------------------------------------------------
  do k=1,nz-1
   do j=js_pe-1,je_pe
    do i=is_pe-1,ie_pe
     fxa = 0.5*(kappaM(i,j,k)+kappaM(i,j+1,k))
     flux_top(i,j,k)=fxa*(v(i,j,k+1,taup1)-v(i,j,k,taup1))/dzw(k)*maskV(i,j,k+1)*maskV(i,j,k)
    enddo
   enddo
  enddo
  do k=1,nz-1
   do j=js_pe-1,je_pe
    do i=is_pe-1,ie_pe
     diss(i,j,k) = (v(i,j,k+1,tau)-v(i,j,k,tau))*flux_top(i,j,k)/dzw(k)  
    enddo
   enddo
  enddo
  diss(:,:,nz)=0.0
  call vgrid_to_tgrid(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,diss,diss)
  K_diss_v = K_diss_v + diss 

 if (.not.enable_hydrostatic) then
  !if (my_pe==0) print'(/a/)','ERROR: implicit vertical friction for vertical velocity not implemented'
  !call halt_stop(' in implicit_vert_friction')

  do j=js_pe,je_pe
   do i=is_pe,ie_pe
    ks=kbot(i,j)
    if (ks>0) then
     do k=ks,nz-1
      delta(k) = dt_mom/dzt(k+1)*0.5*(kappaM(i,j,k)+kappaM(i,j,k+1))
     enddo
     delta(nz)=0.0
     do k=ks+1,nz-1
       a_tri(k) = - delta(k-1)/dzw(k)
     enddo
     a_tri(ks)=0.0
     a_tri(nz) = - delta(nz-1)/(0.5*dzw(nz))
     do k=ks+1,nz-1
      b_tri(k) = 1+ delta(k)/dzw(k) + delta(k-1)/dzw(k) 
     enddo
     b_tri(nz) = 1+ delta(nz-1)/(0.5*dzw(nz))         
     b_tri(ks) = 1+ delta(ks)/dzw(ks)                
     do k=ks,nz-1
      c_tri(k) = - delta(k)/dzw(k)
     enddo
     c_tri(nz)=0.0
     d_tri(ks:nz)=w(i,j,ks:nz,tau) 
     call solve_tridiag(a_tri(ks:nz),b_tri(ks:nz),c_tri(ks:nz),d_tri(ks:nz),w(i,j,ks:nz,taup1),nz-ks+1)
    endif
    dw_mix(i,j,:)=(w(i,j,:,taup1)-w(i,j,:,tau))/dt_mom
   enddo
  enddo

 !---------------------------------------------------------------------------------
 ! diagnose dissipation by vertical friction of vertical momentum
 !---------------------------------------------------------------------------------
  do k=1,nz-1
   do j=js_pe-1,je_pe
    do i=is_pe-1,ie_pe
     fxa = 0.5*(kappaM(i,j,k)+kappaM(i,j,k+1))
     flux_top(i,j,k)=fxa*(w(i,j,k+1,taup1)-w(i,j,k,taup1))/dzt(k+1)*maskW(i,j,k+1)*maskW(i,j,k)
    enddo
   enddo
  enddo
  do k=1,nz-1
   do j=js_pe-1,je_pe
    do i=is_pe-1,ie_pe
     diss(i,j,k) = (w(i,j,k+1,tau)-w(i,j,k,tau))*flux_top(i,j,k)/dzt(k+1)  
    enddo
   enddo
  enddo
  diss(:,:,nz)=0.0
  K_diss_v = K_diss_v + diss 

 endif

end subroutine implicit_vert_friction





