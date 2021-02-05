




subroutine set_idemix3_parameter
!=======================================================================
! set main IDEMIX parameter 
!=======================================================================
  use main_module   
  use idemix_module   
  implicit none
  !include "mass.include"  ! include this on AIX which does not know function acosh, also link with -lmass
  real*8 :: bN0,fxa,Nloc,floc,fxb
  real*8 :: gofx2,hofx1
  integer :: i,j,k
  real*8 :: uz(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz)
  real*8 :: vz(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz)
  real*8 :: aloc(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz)
  real*8 :: yN,yS
  

 !---------------------------------------------------------------------------------
 ! Idemix 1 parameter
 !---------------------------------------------------------------------------------
  do j=js_pe-onx,je_pe+onx
   do i=is_pe-onx,ie_pe+onx
    bN0=0.0
    do k=1,nz-1
     bN0 = bN0 + max(0d0,Nsqr(i,j,k,tau))**0.5*dzw(k)*maskW(i,j,k) 
    enddo
    bN0 = bN0 + max(0d0,Nsqr(i,j,nz,tau))**0.5*0.5*dzw(nz)*maskW(i,j,nz) 
    cstar(i,j) = max(1d-2,bN0/(pi*jstar) )
    do k=1,nz
     fxa = max(0d0,Nsqr(i,j,k,tau))**0.5/(1d-22 + abs(coriolis_t(i,j)) )
     alpha_c(i,j,k) = max( mu0_min, mu0*acosh(max(1d0,fxa))*abs(coriolis_t(i,j))/cstar(i,j)**2 )*maskW(i,j,k) 
     c0(i,j,k)=max(0d0, gamma*cstar(i,j)*gofx2(fxa)*maskW(i,j,k) )
     v0(i,j,k)=max(0d0, gamma*cstar(i,j)*hofx1(fxa)*maskW(i,j,k) )
    enddo
   enddo
  enddo
 !---------------------------------------------------------------------------------
 !---------------------------------------------------------------------------------
  do i=is_pe-onx,ie_pe+onx-1
    do k=1,nz-1
     c0_u(i,:,k) = 0.5*(c0(i,:,k)+c0(i+1,:,k))*maskU(i,:,k)*maskU(i,:,k+1)
     alpha_c_u(i,:,k) = 0.5*(alpha_c(i,:,k)+alpha_c(i+1,:,k))*maskU(i,:,k)*maskU(i,:,k+1)
    enddo
    k=nz
    c0_u(i,:,k) = 0.5*(c0(i,:,k)+c0(i+1,:,k))*maskU(i,:,k)
    alpha_c_u(i,:,k) = 0.5*(alpha_c(i,:,k)+alpha_c(i+1,:,k))*maskU(i,:,k)
  enddo

  do j=js_pe-onx,je_pe+onx-1
    do k=1,nz-1
     c0_v(:,j,k) = 0.5*(c0(:,j,k)+c0(:,j+1,k))*maskV(:,j,k)*maskV(:,j,k+1)
     alpha_c_v(:,j,k) = 0.5*(alpha_c(:,j,k)+alpha_c(:,j+1,k))*maskV(:,j,k)*maskV(:,j,k+1)
    enddo
    k=nz
    c0_v(:,j,k) = 0.5*(c0(:,j,k)+c0(:,j+1,k))*maskV(:,j,k)
    alpha_c_v(:,j,k) = 0.5*(alpha_c(:,j,k)+alpha_c(:,j+1,k))*maskV(:,j,k)
  enddo


 !---------------------------------------------------------------------------------
 !  Shear
 !---------------------------------------------------------------------------------
 do k=1,nz-1
    uz(:,:,k) = (u(:,:,k+1,tau)-u(:,:,k,tau))/dzw(k)*maskU(:,:,k)*maskU(:,:,k+1)
    vz(:,:,k) = (v(:,:,k+1,tau)-v(:,:,k,tau))/dzw(k)*maskV(:,:,k)*maskV(:,:,k+1)
 enddo
 uz(:,:,nz)=uz(:,:,nz-1); vz(:,:,nz)=vz(:,:,nz-1)

 !---------------------------------------------------------------------------------
 ! Idemix 3 parameter
 !---------------------------------------------------------------------------------
  do k=1,nz
   do j=js_pe-onx,je_pe+onx
     do i=is_pe-onx,ie_pe+onx
        floc = max(abs(coriolis_t(i,j)),1D-6)
        Nloc = max(floc*1.2, sqrt( max(0d0, Nsqr(i,j,k,tau))  ) )
        fxa = 0.5*floc/Nloc*(Nloc**2+0.5*floc**2)/(Nloc**2-floc**2 ) 
        fxa = fxa * log( (Nloc+sqrt(Nloc**2-floc**2) )/(Nloc-sqrt(Nloc**2-floc**2)) )
        fxa = fxa - 3./2.*floc/sqrt(Nloc**2-floc**2 ) 
        aloc(i,j,k) = fxa*(2./pi)/(1-2./pi*asin(floc/Nloc) )*maskW(i,j,k)
        fxa=0.5*(1.+tanh( (Nloc/floc-Noverf_min0)/Noverf_min1)  )
        fxb=0.5*(1.+tanh(  (abs(coriolis_t(i,j) ) -1e-6)/1e-8 )     )
        aloc(i,j,k) = fxb*fxa*aloc(i,j,k)
    enddo
   enddo
  enddo

  do j=js_pe-onx,je_pe+onx-1
    do i=is_pe-onx,ie_pe+onx-1
        Om_id3_u(i,j,:) = 0.5*(aloc(i,j,:) + aloc(i+1,j,:) )
        Om_id3_v(i,j,:) = 0.5*(aloc(i,j,:) + aloc(i,j+1,:) )
   enddo
  enddo


 if (enable_idemix3_obc_taper_north) then
    yN = 0.0; if (my_blk_j == n_pes_j) yN = yu(ny)
    call global_max(yN)
    do j=js_pe-onx,je_pe+onx
      Om_id3_u(:,j,:)=  Om_id3_u(:,j,:)*(1- exp( -(yN-yt(j))**2/id3_north_taper_len**2  ) )
      Om_id3_v(:,j,:)=  Om_id3_v(:,j,:)*(1- exp( -(yN-yt(j))**2/id3_north_taper_len**2  ) )
    enddo
 endif
 
 if (enable_idemix3_obc_taper_south) then
    yS = 1e20; if (my_blk_j == 1) yS = yu(0)
    call global_min(yS)
    do j=js_pe-onx,je_pe+onx
      Om_id3_u(:,j,:)=  Om_id3_u(:,j,:)*(1- exp( -(yt(j)-yS)**2/id3_south_taper_len**2  ) )
      Om_id3_v(:,j,:)=  Om_id3_v(:,j,:)*(1- exp( -(yt(j)-yS)**2/id3_south_taper_len**2  ) )
    enddo
 endif


  do k=1,nz
   do j=js_pe-onx,je_pe+onx-1
     do i=is_pe-onx,ie_pe+onx-1
        Om_id3_u(i,j,k) = min(Om_id3_u(i,j,k),  tc_max/sqrt(8.)*pi/(1e-12+abs(uz(i,j,k))) ) ! limit Omega instead of tc_u
        Om_id3_v(i,j,k) = min(Om_id3_v(i,j,k),  tc_max/sqrt(8.)*pi/(1e-12+abs(vz(i,j,k))) )
        Tc_u(i,j,k) = sqrt(8.)/pi*uz(i,j,k)*Om_id3_u(i,j,k) 
        Tc_v(i,j,k) = sqrt(8.)/pi*vz(i,j,k)*Om_id3_v(i,j,k) 
    enddo
   enddo
  enddo

 if (enable_idemix_hor_diffusion_iter) then
    ! check for stability criterium, lateral diffusion is explicit
    !  tau_h v0^2 *dt/dx^2 <= 0.5  ->   v0  <  sqrt( 0.5*dx^2/(dt tau_h)  )
    do j=js_pe-onx,je_pe+onx
      do i=is_pe-onx,ie_pe+onx
         fxa = 0.2*min( dxt(i)*cost(j), dyt(j) )**2/ max(1D0,dt_tracer/idemix_hor_diffusion_iter *tau_h )
         v0(i,j,:) = min( sqrt(fxa), v0(i,j,:) )
      enddo
    enddo
    call set_obc_boundary_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,v0)

    do i=is_pe-onx,ie_pe+onx-1
     v0_u(i,:,:) = 0.5*(v0(i,:,:)+v0(i+1,:,:))*maskU(i,:,:)
    enddo
 
    do j=js_pe-onx,je_pe+onx-1
     v0_v(:,j,:) = 0.5*(v0(:,j,:)+v0(:,j+1,:))*maskV(:,j,:)
    enddo

 endif

end subroutine set_idemix3_parameter



subroutine integrate_idemix3
!=======================================================================
!
!=======================================================================
 use main_module   
 use idemix_module   
 implicit none
 integer :: k,i,j,ks,km2,n
 real*8 :: a_tri(nz),b_tri(nz),c_tri(nz),d_tri(nz),xloc(nz),yloc(nz),zloc(nz)
 real*8 :: forc(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz)
 real*8 :: forc_u(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz)
 real*8 :: forc_v(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz)
 real*8 :: a_w_t(nz),b_w_t(nz),fxa

 call idemix_forcing(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,forc)
 call wgrid_to_tgrid(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,forc,forc)
 call border_exchg_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,forc) 
 call setcyclic_xyz   (is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,forc)
 call tgrid_to_vel_grid(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,forc,forc_u,forc_v)

 !---------------------------------------------------------------------------------
 ! eastward wave propagation F = int_-pi/4^pi/4 ep(phi) dphi
 !
 ! F_t = (tau_v c0 (c0 F)_z)_z  + (tau_v c0 Tc F)_z + tau Tc (c0 F)_z 
 !       + tau Tc^2 F - alpha_c F (F+G)
 !---------------------------------------------------------------------------------
 do j=js_pe,je_pe
  do i=is_pe,ie_pe
    ks=max(kbot(i,j),kbot(i+1,j))
    if (ks>0) then

      a_w_t(ks) = dzw(ks)/(2*dzt(ks))
      b_w_t(ks) = 0
      do k=ks+1,nz-1
        a_w_t(k) = dzw(k)/(2*dzt(k))
        b_w_t(k) = dzw(k-1)/(2*dzt(k))
      enddo  
      a_w_t(nz) = 0.5*dzw(nz)/dzt(nz)
      b_w_t(nz) = dzw(nz-1)/(2*dzt(nz))

      k=ks
      xloc(k) = 1+dt_tracer*alpha_c_u(i,j,k)*max(0d0,G_s(i,j,k,tau)+F_s(i,j,k,tau))
      yloc(k) = c0_u(i,j,k)*dt_tracer/dzt(k) + a_w_t(k)*Tc_u(i,j,k)*dt_tracer
      zloc(k) = c0_u(i,j,k)*dt_tracer/dzt(k) !- b_w_t(k)*Tc_u(i,j,k)*dt_tracer
      do k=ks+1,nz
        xloc(k) = 1+dt_tracer*0.5*(alpha_c_u(i,j,k)+alpha_c_u(i,j,k-1))*max(0d0,G_s(i,j,k,tau)+F_s(i,j,k,tau))
        yloc(k) = c0_u(i,j,k  )*dt_tracer/dzt(k) + a_w_t(k)*Tc_u(i,j,k)*dt_tracer
        zloc(k) = c0_u(i,j,k-1)*dt_tracer/dzt(k) - b_w_t(k)*Tc_u(i,j,k-1)*dt_tracer
      enddo
      k=ks
      a_tri(k) = 0
      b_tri(k) = xloc(k)-tau_v*yloc(k)*Tc_u(i,j,k)*0.5 + tau_v*yloc(k)*c0_u(i,j,k)/(dzw(k))
      c_tri(k) = -tau_v*yloc(k)*(c0_u(i,j,k)+c0_u(i,j,k+1))/(2*dzw(k)) - tau_v*yloc(k)*Tc_u(i,j,k)*0.5
      do k=ks+1,nz-1
        km2=max(ks,k-2)
        a_tri(k) = -tau_v*zloc(k)*(c0_u(i,j,km2)+c0_u(i,j,k-1))/(2*dzw(k-1)) + tau_v*zloc(k)*Tc_u(i,j,k-1)*0.5
        b_tri(k) = xloc(k) - tau_v*yloc(k)*Tc_u(i,j,k)*0.5 + tau_v*zloc(k)*Tc_u(i,j,k-1)*0.5 &
                           + tau_v*yloc(k)*(c0_u(i,j,k-1)+c0_u(i,j,k))/(2*dzw(k))  &
                           + tau_v*zloc(k)*(c0_u(i,j,k-1)+c0_u(i,j,k))/(2*dzw(k-1))
        c_tri(k) = -tau_v*yloc(k)*(c0_u(i,j,k)+c0_u(i,j,k+1))/(2*dzw(k)) - tau_v*yloc(k)*Tc_u(i,j,k)*0.5
      enddo
      k=nz
      a_tri(k) = -tau_v*zloc(k)*(c0_u(i,j,k-2)+c0_u(i,j,k-1))/(2*dzw(k-1)) + tau_v*zloc(k)*Tc_u(i,j,k-1)*0.5
      b_tri(k) = xloc(k)+tau_v*zloc(k)*Tc_u(i,j,k-1)*0.5 + tau_v*zloc(k)*(c0_u(i,j,k-1)+c0_u(i,j,k))/(2*dzw(k-1))
      c_tri(k) = 0

      d_tri(ks:nz) = F_s(i,j,ks:nz,tau) + dt_tracer*(0.5*forc_u(i,j,ks:nz)  )
      !d_tri(ks) = d_tri(ks) + zloc(ks)/c0_u(i,j,ks)*0.25*forc_iw_bottom_u(i,j)
      !d_tri(nz) = d_tri(nz) + yloc(nz)/c0_u(i,j,nz)*0.25*forc_iw_surface_u(i,j)
      d_tri(ks) = d_tri(ks) + dt_tracer/dzt(ks)*0.25*forc_iw_bottom_u(i,j)    
      d_tri(nz) = d_tri(nz) + dt_tracer/dzt(nz)*0.25*forc_iw_surface_u(i,j)
      call solve_tridiag(a_tri(ks:nz),b_tri(ks:nz),c_tri(ks:nz),d_tri(ks:nz),F_s(i,j,ks:nz,taup1),nz-ks+1)
    endif
  enddo
 enddo


 !---------------------------------------------------------------------------------
 ! westward wave propagation G = int_3pi/4^5pi/4 ep(phi) dphi
 !
 ! G_t = (tau_v c0 (c0 G)_z)_z  - (tau_v c0 Tc G)_z - tau Tc (c0 G)_z 
 !       + tau Tc^2 G - alpha_c G (F+G)
 !---------------------------------------------------------------------------------
 do j=js_pe,je_pe
  do i=is_pe,ie_pe
    ks=max(kbot(i,j),kbot(i+1,j))
    if (ks>0) then

      a_w_t(ks) = dzw(ks)/(2*dzt(ks))
      b_w_t(ks) = 0
      do k=ks+1,nz-1
        a_w_t(k) = dzw(k)/(2*dzt(k))
        b_w_t(k) = dzw(k-1)/(2*dzt(k))
      enddo  
      a_w_t(nz) = 0.5*dzw(nz)/dzt(nz)
      b_w_t(nz) = dzw(nz-1)/(2*dzt(nz))

      k=ks
      xloc(k) = 1+dt_tracer*alpha_c_u(i,j,k)*max(0d0,G_s(i,j,k,tau)+F_s(i,j,k,tau))
      yloc(k) = c0_u(i,j,k)*dt_tracer/dzt(k) - a_w_t(k)*Tc_u(i,j,k)*dt_tracer
      zloc(k) = c0_u(i,j,k)*dt_tracer/dzt(k) !+ b_w_t(k)*Tc_u(i,j,k)*dt_tracer
      do k=ks+1,nz
        xloc(k) = 1+0.5*dt_tracer*(alpha_c_u(i,j,k)+alpha_c_u(i,j,k-1))*max(0d0,G_s(i,j,k,tau)+F_s(i,j,k,tau))
        yloc(k) = c0_u(i,j,k  )*dt_tracer/dzt(k) - a_w_t(k)*Tc_u(i,j,k)*dt_tracer
        zloc(k) = c0_u(i,j,k-1)*dt_tracer/dzt(k) + b_w_t(k)*Tc_u(i,j,k-1)*dt_tracer
      enddo
      k=ks
      a_tri(k) = 0
      b_tri(k) = xloc(k)+tau_v*yloc(k)*Tc_u(i,j,k)*0.5 + tau_v*yloc(k)*c0_u(i,j,k)/(dzw(k))
      c_tri(k) = -tau_v*yloc(k)*(c0_u(i,j,k)+c0_u(i,j,k+1))/(2*dzw(k)) + tau_v*yloc(k)*Tc_u(i,j,k)*0.5
      do k=ks+1,nz-1
        km2=max(ks,k-2)
        a_tri(k) = -tau_v*zloc(k)*(c0_u(i,j,km2)+c0_u(i,j,k-1))/(2*dzw(k-1)) - tau_v*zloc(k)*Tc_u(i,j,k-1)*0.5
        b_tri(k) = xloc(k) + tau_v*yloc(k)*Tc_u(i,j,k)*0.5 - tau_v*zloc(k)*Tc_u(i,j,k-1)*0.5 &
                           + tau_v*yloc(k)*(c0_u(i,j,k-1)+c0_u(i,j,k))/(2*dzw(k))  &
                           + tau_v*zloc(k)*(c0_u(i,j,k-1)+c0_u(i,j,k))/(2*dzw(k-1))
        c_tri(k) = -tau_v*yloc(k)*(c0_u(i,j,k)+c0_u(i,j,k+1))/(2*dzw(k)) + tau_v*yloc(k)*Tc_u(i,j,k)*0.5
      enddo
      k=nz
      a_tri(k) = -tau_v*zloc(k)*(c0_u(i,j,k-2)+c0_u(i,j,k-1))/(2*dzw(k-1)) - tau_v*zloc(k)*Tc_u(i,j,k-1)*0.5
      b_tri(k) = xloc(k)-tau_v*zloc(k)*Tc_u(i,j,k-1)*0.5 + tau_v*zloc(k)*(c0_u(i,j,k-1)+c0_u(i,j,k))/(2*dzw(k-1))
      c_tri(k) = 0

      d_tri(ks:nz) = G_s(i,j,ks:nz,tau)+ dt_tracer*(0.5*forc_u(i,j,ks:nz) )
      !d_tri(ks) = d_tri(ks) + zloc(ks)/c0_u(i,j,ks)*0.25*forc_iw_bottom_u(i,j)
      !d_tri(nz) = d_tri(nz) + yloc(nz)/c0_u(i,j,nz)*0.25*forc_iw_surface_u(i,j)
       
      d_tri(ks) = d_tri(ks) + dt_tracer/dzt(ks)*0.25*forc_iw_bottom_u(i,j)
      d_tri(nz) = d_tri(nz) + dt_tracer/dzt(nz)*0.25*forc_iw_surface_u(i,j)
      call solve_tridiag(a_tri(ks:nz),b_tri(ks:nz),c_tri(ks:nz),d_tri(ks:nz),G_s(i,j,ks:nz,taup1),nz-ks+1)
    endif
  enddo
 enddo


 !---------------------------------------------------------------------------------
 ! northward wave propagation O = int_pi/4^3pi/4 ep(phi) dphi
 !---------------------------------------------------------------------------------
 do j=js_pe,je_pe
  do i=is_pe,ie_pe
    ks=max(kbot(i,j),kbot(i,j+1))
    if (ks>0) then

      a_w_t(ks) = dzw(ks)/(2*dzt(ks))
      b_w_t(ks) = 0
      do k=ks+1,nz-1
        a_w_t(k) = dzw(k)/(2*dzt(k))
        b_w_t(k) = dzw(k-1)/(2*dzt(k))
      enddo  
      a_w_t(nz) = 0.5*dzw(nz)/dzt(nz)
      b_w_t(nz) = dzw(nz-1)/(2*dzt(nz))

      k=ks
      xloc(k) = 1+dt_tracer*alpha_c_v(i,j,k)*max(0d0,O_s(i,j,k,tau)+P_s(i,j,k,tau))
      yloc(k) = c0_v(i,j,k)*dt_tracer/dzt(k) + a_w_t(k)*Tc_v(i,j,k)*dt_tracer
      zloc(k) = c0_v(i,j,k)*dt_tracer/dzt(k) !- b_w_t(k)*Tc_v(i,j,k)*dt_tracer
      do k=ks+1,nz
        xloc(k) = 1+dt_tracer*0.5*(alpha_c_v(i,j,k)+alpha_c_v(i,j,k-1))*max(0d0,O_s(i,j,k,tau)+P_s(i,j,k,tau))
        yloc(k) = c0_v(i,j,k  )*dt_tracer/dzt(k) + a_w_t(k)*Tc_v(i,j,k)*dt_tracer
        zloc(k) = c0_v(i,j,k-1)*dt_tracer/dzt(k) - b_w_t(k)*Tc_v(i,j,k-1)*dt_tracer
      enddo
      k=ks
      a_tri(k) = 0
      b_tri(k) = xloc(k)-tau_v*yloc(k)*Tc_v(i,j,k)*0.5 + tau_v*yloc(k)*c0_v(i,j,k)/(dzw(k))
      c_tri(k) = -tau_v*yloc(k)*(c0_v(i,j,k)+c0_v(i,j,k+1))/(2*dzw(k)) - tau_v*yloc(k)*Tc_v(i,j,k)*0.5
      do k=ks+1,nz-1
        km2=max(ks,k-2)
        a_tri(k) = -tau_v*zloc(k)*(c0_v(i,j,km2)+c0_v(i,j,k-1))/(2*dzw(k-1)) + tau_v*zloc(k)*Tc_v(i,j,k-1)*0.5
        b_tri(k) = xloc(k) - tau_v*yloc(k)*Tc_v(i,j,k)*0.5 + tau_v*zloc(k)*Tc_v(i,j,k-1)*0.5 &
                           + tau_v*yloc(k)*(c0_v(i,j,k-1)+c0_v(i,j,k))/(2*dzw(k))  &
                           + tau_v*zloc(k)*(c0_v(i,j,k-1)+c0_v(i,j,k))/(2*dzw(k-1))
        c_tri(k) = -tau_v*yloc(k)*(c0_v(i,j,k)+c0_v(i,j,k+1))/(2*dzw(k)) - tau_v*yloc(k)*Tc_v(i,j,k)*0.5
      enddo
      k=nz
      a_tri(k) = -tau_v*zloc(k)*(c0_v(i,j,k-2)+c0_v(i,j,k-1))/(2*dzw(k-1)) + tau_v*zloc(k)*Tc_v(i,j,k-1)*0.5
      b_tri(k) = xloc(k)+tau_v*zloc(k)*Tc_v(i,j,k-1)*0.5 + tau_v*zloc(k)*(c0_v(i,j,k-1)+c0_v(i,j,k))/(2*dzw(k-1))
      c_tri(k) = 0

      d_tri(ks:nz) = O_s(i,j,ks:nz,tau) + dt_tracer*(0.5*forc_v(i,j,ks:nz)  )
      !d_tri(ks) = d_tri(ks) + zloc(ks)/c0_v(i,j,ks)*0.25*forc_iw_bottom_v(i,j) 
      !d_tri(nz) = d_tri(nz) + yloc(nz)/c0_v(i,j,nz)*0.25*forc_iw_surface_v(i,j)
      d_tri(ks) = d_tri(ks) + dt_tracer/dzt(ks)*0.25*forc_iw_bottom_v(i,j)
      d_tri(nz) = d_tri(nz) + dt_tracer/dzt(nz)*0.25*forc_iw_surface_v(i,j)
      call solve_tridiag(a_tri(ks:nz),b_tri(ks:nz),c_tri(ks:nz),d_tri(ks:nz),O_s(i,j,ks:nz,taup1),nz-ks+1)
    endif
  enddo
 enddo

 !---------------------------------------------------------------------------------
 ! southward wave propagation P = int_5pi/4^7pi/4 ep(phi) dphi
 !---------------------------------------------------------------------------------
 do j=js_pe,je_pe
  do i=is_pe,ie_pe
    ks=max(kbot(i,j),kbot(i,j+1))
    if (ks>0) then

      a_w_t(ks) = dzw(ks)/(2*dzt(ks))
      b_w_t(ks) = 0
      do k=ks+1,nz-1
        a_w_t(k) = dzw(k)/(2*dzt(k))
        b_w_t(k) = dzw(k-1)/(2*dzt(k))
      enddo  
      a_w_t(nz) = 0.5*dzw(nz)/dzt(nz)
      b_w_t(nz) = dzw(nz-1)/(2*dzt(nz))

      k=ks
      xloc(k) = 1+dt_tracer*alpha_c_v(i,j,k)*max(0d0,P_s(i,j,k,tau)+O_s(i,j,k,tau))
      yloc(k) = c0_v(i,j,k)*dt_tracer/dzt(k) - a_w_t(k)*Tc_v(i,j,k)*dt_tracer
      zloc(k) = c0_v(i,j,k)*dt_tracer/dzt(k) !+ b_w_t(k)*Tc_v(i,j,k)*dt_tracer
      do k=ks+1,nz
        xloc(k) = 1+0.5*dt_tracer*(alpha_c_v(i,j,k)+alpha_c_v(i,j,k-1))*max(0d0,O_s(i,j,k,tau)+P_s(i,j,k,tau))
        yloc(k) = c0_v(i,j,k  )*dt_tracer/dzt(k) - a_w_t(k)*Tc_v(i,j,k)*dt_tracer
        zloc(k) = c0_v(i,j,k-1)*dt_tracer/dzt(k) + b_w_t(k)*Tc_v(i,j,k-1)*dt_tracer
      enddo
      k=ks
      a_tri(k) = 0
      b_tri(k) = xloc(k)+tau_v*yloc(k)*Tc_v(i,j,k)*0.5 + tau_v*yloc(k)*c0_v(i,j,k)/(dzw(k))
      c_tri(k) = -tau_v*yloc(k)*(c0_v(i,j,k)+c0_v(i,j,k+1))/(2*dzw(k)) + tau_v*yloc(k)*Tc_v(i,j,k)*0.5
      do k=ks+1,nz-1
        km2=max(ks,k-2)
        a_tri(k) = -tau_v*zloc(k)*(c0_v(i,j,km2)+c0_v(i,j,k-1))/(2*dzw(k-1)) - tau_v*zloc(k)*Tc_v(i,j,k-1)*0.5
        b_tri(k) = xloc(k) + tau_v*yloc(k)*Tc_v(i,j,k)*0.5 - tau_v*zloc(k)*Tc_v(i,j,k-1)*0.5 &
                           + tau_v*yloc(k)*(c0_v(i,j,k-1)+c0_v(i,j,k))/(2*dzw(k))  &
                           + tau_v*zloc(k)*(c0_v(i,j,k-1)+c0_v(i,j,k))/(2*dzw(k-1))
        c_tri(k) = -tau_v*yloc(k)*(c0_v(i,j,k)+c0_v(i,j,k+1))/(2*dzw(k)) + tau_v*yloc(k)*Tc_v(i,j,k)*0.5
      enddo
      k=nz
      a_tri(k) = -tau_v*zloc(k)*(c0_v(i,j,k-2)+c0_v(i,j,k-1))/(2*dzw(k-1)) - tau_v*zloc(k)*Tc_v(i,j,k-1)*0.5
      b_tri(k) = xloc(k)-tau_v*zloc(k)*Tc_v(i,j,k-1)*0.5 + tau_v*zloc(k)*(c0_v(i,j,k-1)+c0_v(i,j,k))/(2*dzw(k-1))
      c_tri(k) = 0

      d_tri(ks:nz) = P_s(i,j,ks:nz,tau)+ dt_tracer*(0.5*forc_v(i,j,ks:nz) )
      !d_tri(ks) = d_tri(ks) + zloc(ks)/c0_v(i,j,ks)*0.25*forc_iw_bottom_v(i,j)
      !d_tri(nz) = d_tri(nz) + yloc(nz)/c0_v(i,j,nz)*0.25*forc_iw_surface_v(i,j)
      d_tri(ks) = d_tri(ks) + dt_tracer/dzt(ks)*0.25*forc_iw_bottom_v(i,j)
      d_tri(nz) = d_tri(nz) + dt_tracer/dzt(nz)*0.25*forc_iw_surface_v(i,j)
      call solve_tridiag(a_tri(ks:nz),b_tri(ks:nz),c_tri(ks:nz),d_tri(ks:nz),P_s(i,j,ks:nz,taup1),nz-ks+1)
    endif
  enddo
 enddo

 if (enable_idemix3_cut_negative_energy) then
   ! very dangerous switch. Use with care
   F_s(is_pe:ie_pe,js_pe:je_pe,:,taup1) = max(0d0,F_s(is_pe:ie_pe,js_pe:je_pe,:,taup1))
   G_s(is_pe:ie_pe,js_pe:je_pe,:,taup1) = max(0d0,G_s(is_pe:ie_pe,js_pe:je_pe,:,taup1))
   O_s(is_pe:ie_pe,js_pe:je_pe,:,taup1) = max(0d0,O_s(is_pe:ie_pe,js_pe:je_pe,:,taup1))
   P_s(is_pe:ie_pe,js_pe:je_pe,:,taup1) = max(0d0,P_s(is_pe:ie_pe,js_pe:je_pe,:,taup1))
 endif

 call border_exchg_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,F_s(:,:,:,taup1)) 
 call setcyclic_xyz   (is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,F_s(:,:,:,taup1))
 call set_obc_boundary_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,F_s(:,:,:,taup1))

 call border_exchg_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,G_s(:,:,:,taup1)) 
 call setcyclic_xyz   (is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,G_s(:,:,:,taup1))
 call set_obc_boundary_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,G_s(:,:,:,taup1))

 call border_exchg_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,O_s(:,:,:,taup1)) 
 call setcyclic_xyz   (is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,O_s(:,:,:,taup1))
 call set_obc_boundary_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,O_s(:,:,:,taup1))

 call border_exchg_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,P_s(:,:,:,taup1)) 
 call setcyclic_xyz   (is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,P_s(:,:,:,taup1))
 call set_obc_boundary_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,P_s(:,:,:,taup1))

 !---------------------------------------------------------------------------------
 ! Diagnose Dissipation 
 !---------------------------------------------------------------------------------
 forc_u = 0.; forc_v = 0.
 do j=js_pe-onx,je_pe+onx-1
  do i=is_pe-onx,ie_pe+onx-1
    ks=max(kbot(i,j),kbot(i+1,j))
    if (ks>0) then
      k=ks
      forc_u(i,j,k) = alpha_c_u(i,j,k)*max(0d0,G_s(i,j,k,tau)+F_s(i,j,k,tau))*(G_s(i,j,k,taup1)+F_s(i,j,k,taup1)) *maskU(i,j,k)
      do k=ks+1,nz
        forc_u(i,j,k) = 0.5*(alpha_c_u(i,j,k)+alpha_c_u(i,j,k-1))*  &
                         max(0d0,G_s(i,j,k,tau)+F_s(i,j,k,tau))*(G_s(i,j,k,taup1)+F_s(i,j,k,taup1)) *maskU(i,j,k)
      enddo
    endif
    ks=max(kbot(i,j),kbot(i,j+1))
    if (ks>0) then
      k=ks
      forc_v(i,j,k) = alpha_c_v(i,j,k)*max(0d0,O_s(i,j,k,tau)+P_s(i,j,k,tau))*(O_s(i,j,k,taup1)+P_s(i,j,k,taup1)) *maskV(i,j,k)
      do k=ks+1,nz
        forc_v(i,j,k) = 0.5*(alpha_c_v(i,j,k)+alpha_c_v(i,j,k-1))*  &
                         max(0d0,O_s(i,j,k,tau)+P_s(i,j,k,tau))*(O_s(i,j,k,taup1)+P_s(i,j,k,taup1)) *maskV(i,j,k)
      enddo
    endif
  enddo
 enddo
 call ugrid_to_tgrid(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,forc_u,forc_u)
 call vgrid_to_tgrid(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,forc_v,forc_v)
 iw_diss = forc_u +forc_v
 call tgrid_to_wgrid(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,iw_diss,iw_diss)


 if (enable_idemix_hor_diffusion_iter) then
  !---------------------------------------------------------------------------------
  ! add tendency due to lateral diffusion with iterative method
  !---------------------------------------------------------------------------------

  fxa = dt_tracer/idemix_hor_diffusion_iter
  do n=1,idemix_hor_diffusion_iter

   call idemix3_diffusion_u(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,F_s(:,:,:,taup1),fxa)
   call border_exchg_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,F_s(:,:,:,taup1)) 
   call setcyclic_xyz   (is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,F_s(:,:,:,taup1))
   call set_obc_boundary_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,F_s(:,:,:,taup1))

   call idemix3_diffusion_u(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,G_s(:,:,:,taup1),fxa)
   call border_exchg_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,G_s(:,:,:,taup1)) 
   call setcyclic_xyz   (is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,G_s(:,:,:,taup1))
   call set_obc_boundary_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,G_s(:,:,:,taup1))

   call idemix3_diffusion_v(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,O_s(:,:,:,taup1),fxa)
   call border_exchg_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,O_s(:,:,:,taup1)) 
   call setcyclic_xyz   (is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,O_s(:,:,:,taup1))
   call set_obc_boundary_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,O_s(:,:,:,taup1))

   call idemix3_diffusion_v(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,P_s(:,:,:,taup1),fxa)
   call border_exchg_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,P_s(:,:,:,taup1)) 
   call setcyclic_xyz   (is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,P_s(:,:,:,taup1))
   call set_obc_boundary_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,P_s(:,:,:,taup1))

  enddo

 endif

 call idemix_get_delta(taup1)
 !call idemix_friction(taup1)

end subroutine integrate_idemix3






subroutine idemix3_diffusion_u(is_,ie_,js_,je_,nz_,Aloc,dtloc)
!--------------------------------------------------------------
!--------------------------------------------------------------
 use main_module   
 use idemix_module   
 implicit none
 integer, intent(in) :: is_,ie_,js_,je_,nz_
 real*8, intent(inout)  :: Aloc(is_:ie_,js_:je_,nz_)
 real*8, intent(in) :: dtloc
 integer :: i,j

 do j=js_pe,je_pe
   do i=is_pe-1,ie_pe
    flux_east(i,j,:)=tau_h*0.5*(v0_u(i+1,j,:)+v0_u(i,j,:))*( v0_u(i+1,j,:)*Aloc(i+1,j,:) &
                                -v0_u(i,j,:)*Aloc(i,j,:))/(cost(j)*dxt(i+1))*maskU(i+1,j,:)*maskU(i,j,:)
   enddo
 enddo
 do j=js_pe-1,je_pe
    flux_north(:,j,:)=tau_h*0.5*(v0_u(:,j,:)+v0_u(:,j+1,:))*( v0_u(:,j+1,:)*Aloc(:,j+1,:) &
                                         -v0_u(:,j,:)*Aloc(:,j,:))/dyu(j)*maskU(:,j+1,:)*maskU(:,j,:)*cosu(j)
 enddo 

 !---------------------------------------------------------------------------------
 ! update tendency 
 !---------------------------------------------------------------------------------
 do j=js_pe,je_pe
   do i=is_pe,ie_pe
    Aloc(i,j,:)= Aloc(i,j,:) + dtloc*maskU(i,j,:)*((flux_east(i,j,:)  - flux_east(i-1,j,:))/(cost(j)*dxu(i)) &
                                                  +(flux_north(i,j,:) - flux_north(i,j-1,:))/(cost(j)*dyt(j)) )
  enddo
 enddo
end subroutine idemix3_diffusion_u





subroutine idemix3_diffusion_v(is_,ie_,js_,je_,nz_,Aloc,dtloc)
!--------------------------------------------------------------
!--------------------------------------------------------------
 use main_module   
 use idemix_module   
 implicit none
 integer, intent(in) :: is_,ie_,js_,je_,nz_
 real*8, intent(inout)  :: Aloc(is_:ie_,js_:je_,nz_)
 real*8, intent(in) :: dtloc
 integer :: i,j

 do j=js_pe,je_pe
   do i=is_pe-1,ie_pe
    flux_east(i,j,:)=tau_h*0.5*(v0_v(i,j,:)+v0_v(i+1,j,:))*( v0_v(i+1,j,:)*Aloc(i+1,j,:) &
                                        -v0_v(i,j,:)*Aloc(i,j,:))/(cosu(j)*dxu(i))*maskV(i+1,j,:)*maskV(i,j,:)
   enddo
 enddo
 do j=js_pe-1,je_pe
    flux_north(:,j,:)=tau_h*0.5*(v0_v(:,j,:)+v0_v(:,j+1,:))*( v0_v(:,j+1,:)*Aloc(:,j+1,:) &
                                        -v0_v(:,j,:)*Aloc(:,j,:))/dyt(j+1)*maskV(:,j+1,:)*maskV(:,j,:)*cost(j+1)
 enddo 

 !---------------------------------------------------------------------------------
 ! update tendency 
 !---------------------------------------------------------------------------------
 do j=js_pe,je_pe
   do i=is_pe,ie_pe
    Aloc(i,j,:)= Aloc(i,j,:) + dtloc*maskV(i,j,:)*((flux_east(i,j,:)  - flux_east(i-1,j,:))/(cosu(j)*dxt(i)) &
                                                  +(flux_north(i,j,:) - flux_north(i,j-1,:))/(cosu(j)*dyu(j)) )
  enddo
 enddo
end subroutine idemix3_diffusion_v
