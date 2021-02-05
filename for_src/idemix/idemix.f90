


subroutine set_idemix_parameter
  use main_module   
  use idemix_module   
  if (enable_idemix) then
     if (enable_idemix3) then
        call set_idemix3_parameter
     else
        call set_idemix1_parameter
      endif
      if (enable_idemix_M2 .or. enable_idemix_niw) call set_spectral_parameter
  endif
end subroutine set_idemix_parameter


subroutine integrate_idemix
  use main_module   
  use idemix_module   
  use timing_module   
  
  if (enable_leewaves) then

     call tic('leewave_flux')
     call leewave_flux 
     call toc('leewave_flux') 
     
     call tic('integrate_leewaves')
     call integrate_leewaves
     call toc('integrate_leewaves')
  endif
  
  if (enable_idemix_M2)   call integrate_idemix_M2
  if (enable_idemix_niw)  call integrate_idemix_niw
  if (enable_idemix)  then
       if (enable_idemix3)  then
            if (.not. enable_momentum_equation) call idemix_friction
            call integrate_idemix3
       else
            call integrate_idemix1
       endif
  endif
  if (enable_idemix_M2 .or. enable_idemix_niw ) call wave_interaction
end subroutine integrate_idemix


subroutine set_idemix1_parameter
!=======================================================================
! set IDEMIX 1.0 parameter 
!=======================================================================
  use main_module   
  use idemix_module   
  implicit none
  integer :: i,j,k
  real*8 :: fxa,gofx2,bN0,hofx1
  real*8 :: aloc(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz)

  !include "mass.include"  ! include this on AIX which does not know function
  do k=1,nz
   do j=js_pe-onx,je_pe+onx
    do i=is_pe-onx,ie_pe+onx
      aloc(i,j,k) = max(0d0,Nsqr(i,j,k,tau))**0.5
     enddo
   enddo
  enddo
  
  do j=js_pe-onx,je_pe+onx
   do i=is_pe-onx,ie_pe+onx
    bN0=0.0
    do k=1,nz-1
     bN0 = bN0 + aloc(i,j,k)*dzw(k)*maskW(i,j,k) 
    enddo
    bN0 = bN0 + aloc(i,j,nz)*0.5*dzw(nz)*maskW(i,j,nz) 
    do k=1,nz
     fxa = aloc(i,j,k)/(1d-22 + abs(coriolis_t(i,j)) )
     cstar(i,j) = max(1d-2,bN0/(pi*jstar) )
     c0(i,j,k)=max(0d0, gamma*cstar(i,j)*gofx2(fxa)*maskW(i,j,k) )
     v0(i,j,k)=max(0d0, gamma*cstar(i,j)*hofx1(fxa)*maskW(i,j,k) )
     alpha_c(i,j,k) = max( 1d-4,mu0*acosh(max(1d0,fxa))*abs(coriolis_t(i,j))/cstar(i,j)**2 )*maskW(i,j,k) 
    enddo
   enddo
  enddo

  if (enable_idemix_hor_diffusion .or. enable_idemix_hor_diffusion_iter) then
    ! check for stability criterium, lateral diffusion is explicit
    !  tau_h v0^2 *dt/dx^2 <= 0.5  ->   v0  <  sqrt( 0.5*dx^2/(dt tau_h)  )
    do j=js_pe-onx,je_pe+onx
      do i=is_pe-onx,ie_pe+onx
         if (enable_idemix_hor_diffusion_iter) then
            fxa = 0.2*min( dxt(i)*cost(j), dyt(j) )**2/ max(1D0,dt_tracer/idemix_hor_diffusion_iter *tau_h )
         else
            fxa = 0.2*min( dxt(i)*cost(j), dyt(j) )**2/ max(1D0,dt_tracer *tau_h )
         endif
         v0(i,j,:) = min( sqrt(fxa), v0(i,j,:) )
      enddo
    enddo
    call set_obc_boundary_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,v0)
  endif

end subroutine set_idemix1_parameter




subroutine idemix_forcing(is_,ie_,js_,je_,nz_,forc)
 use main_module   
 use eke_module   
 use idemix_module   
 implicit none
 integer :: i,j,k,is_,ie_,js_,je_,nz_,ks
 real*8 :: forc(is_:ie_,js_:je_,nz_)
 real*8 :: a_loc(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx)
 !real*8 :: b_loc(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz)

 !---------------------------------------------------------------------------------
 ! forcing by EKE dissipation
 !---------------------------------------------------------------------------------
 if (enable_eke) then
      forc = eke_diss_iw
 else ! short cut without EKE model
      forc = K_diss_gm  - P_diss_skew
      if (.not. enable_store_lateral_friction_heat)  forc = forc + K_diss_h
      if (.not. enable_store_cabbeling_heat)         forc = forc - P_diss_hmix  - P_diss_iso
 endif

 if (enable_eke.and. (enable_eke_diss_bottom.or.enable_eke_diss_surfbot)) then
 !---------------------------------------------------------------------------------
 ! vertically integrate EKE dissipation and inject at bottom and/or surface
 !---------------------------------------------------------------------------------
   a_loc = 0d0
   do k=1,nz-1
     a_loc = a_loc + dzw(k)*forc(:,:,k)*maskW(:,:,k)
   enddo
   k=nz; a_loc = a_loc + 0.5*dzw(k)*forc(:,:,k)*maskW(:,:,k)
   forc = 0 
   if (enable_eke_diss_bottom) then 
    do j=js_pe,je_pe
     do i=is_pe,ie_pe
      ks = max(1,kbot(i,j))
      forc(i,j,ks)=a_loc(i,j)/dzw(ks)
     enddo
    enddo
   else
    do j=js_pe,je_pe
     do i=is_pe,ie_pe
      ks = max(1,kbot(i,j))
      forc(i,j,ks) =      eke_diss_surfbot_frac*a_loc(i,j)/dzw(ks)
      forc(i,j,nz) = (1.-eke_diss_surfbot_frac)*a_loc(i,j)/(0.5*dzw(nz))
     enddo
    enddo
   endif
 endif

 !---------------------------------------------------------------------------------
 ! forcing by bottom friction
 !---------------------------------------------------------------------------------
 if (.not.enable_store_bottom_friction_tke) forc = forc + K_diss_bot
 
 
 !---------------------------------------------------------------------------------
 ! Lee wave forcing
 !---------------------------------------------------------------------------------
 if (enable_leewaves) then
   if (enable_idemix3) then
     if (my_pe == 0 ) print*,' Error: leewaves with idemix3 not yet implemented'
     call halt_stop('in idemix_forcing')
   endif
   forc = forc + iw_diss_lee
 endif
end subroutine idemix_forcing




subroutine integrate_idemix1
!=======================================================================
! integrate idemix on W grid
!=======================================================================
 use main_module   
 use idemix_module   
 implicit none
 integer :: i,j,k,ks,ke,n
 real*8 :: a_tri(nz),b_tri(nz),c_tri(nz),d_tri(nz),delta(nz)
 real*8 :: forc(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz)
 real*8 :: maxE_iw(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz),fxa
 !real*8 :: a_loc(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx)
 real*8 :: maskUtr,maskVtr
 maskUtr(i,j,k) = maskW(i+1,j,k)*maskW(i,j,k)
 maskVtr(i,j,k) = maskW(i,j+1,k)*maskW(i,j,k)
 
 ke = nz

 call idemix_forcing(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,forc)

 !---------------------------------------------------------------------------------
 ! forcing by EKE dissipation
 !---------------------------------------------------------------------------------
 !if (enable_eke) then
 !     forc = eke_diss_iw
 !else ! short cut without EKE model
 !     forc = K_diss_h + K_diss_gm  - P_diss_skew
 !     if (.not. enable_store_cabbeling_heat)  forc = forc - P_diss_hmix  - P_diss_iso
 !endif

 !if (enable_eke.and. (enable_eke_diss_bottom.or.enable_eke_diss_surfbot)) then
 !---------------------------------------------------------------------------------
 ! vertically integrate EKE dissipation and inject at bottom and/or surface
 !---------------------------------------------------------------------------------
 !  a_loc = 0d0
 !  do k=1,nz-1
 !    a_loc = a_loc + dzw(k)*forc(:,:,k)*maskW(:,:,k)
 !  enddo
 !  k=nz; a_loc = a_loc + 0.5*dzw(k)*forc(:,:,k)*maskW(:,:,k)
 !  forc = 0 
 !  if (enable_eke_diss_bottom) then 
 !   do j=js_pe,je_pe
 !    do i=is_pe,ie_pe
 !     ks = max(1,kbot(i,j))
 !     forc(i,j,ks)=a_loc(i,j)/dzw(ks)
 !    enddo
 !   enddo
 !  else
 !   do j=js_pe,je_pe
 !    do i=is_pe,ie_pe
 !     ks = max(1,kbot(i,j))
 !     forc(i,j,ks) =      eke_diss_surfbot_frac*a_loc(i,j)/dzw(ks)
 !     forc(i,j,ke) = (1.-eke_diss_surfbot_frac)*a_loc(i,j)/(0.5*dzw(ke))
 !    enddo
 !   enddo
 !  endif
 !endif

 !---------------------------------------------------------------------------------
 ! forcing by bottom friction
 !---------------------------------------------------------------------------------
 !if (.not.enable_store_bottom_friction_tke) forc = forc + K_diss_bot



 !---------------------------------------------------------------------------------
 !prevent negative dissipation of IW energy
 !---------------------------------------------------------------------------------
 maxE_iw = max(0D0, E_iw(:,:,:,tau) )

 !---------------------------------------------------------------------------------
 ! vertical diffusion and dissipation is solved implicitely
 !---------------------------------------------------------------------------------
 do j=js_pe,je_pe
  do i=is_pe,ie_pe
    ks=kbot(i,j)
    if (ks>0) then
     do k=ks,ke-1
      delta(k) = dt_tracer*tau_v/dzt(k+1)*0.5*(c0(i,j,k)+c0(i,j,k+1))
     enddo
     delta(ke)=0.0
     do k=ks+1,ke-1
       a_tri(k) = - delta(k-1)*c0(i,j,k-1)/dzw(k)
     enddo
     a_tri(ks)=0.0
     a_tri(ke) = - delta(ke-1)/(0.5*dzw(ke))*c0(i,j,ke-1)
     do k=ks+1,ke-1
      b_tri(k) = 1+ delta(k)*c0(i,j,k)/dzw(k) + delta(k-1)*c0(i,j,k)/dzw(k)+ dt_tracer*alpha_c(i,j,k)*maxE_iw(i,j,k)
     enddo
     b_tri(ke) = 1+ delta(ke-1)/(0.5*dzw(ke))*c0(i,j,ke) + dt_tracer*alpha_c(i,j,ke)*maxE_iw(i,j,ke)
     b_tri(ks) = 1+ delta(ks)/dzw(ks)*c0(i,j,ks)         + dt_tracer*alpha_c(i,j,ks)*maxE_iw(i,j,ks)
     do k=ks,ke-1
      c_tri(k) = - delta(k)/dzw(k)*c0(i,j,k+1)
     enddo
     c_tri(ke) = 0.0
     d_tri(ks:ke) = E_iw(i,j,ks:ke,tau)  + dt_tracer*forc(i,j,ks:ke)
     d_tri(ks) = d_tri(ks) + dt_tracer*forc_iw_bottom(i,j)/dzw(ks)    
     d_tri(ke) = d_tri(ke) + dt_tracer*forc_iw_surface(i,j)/(0.5*dzw(ke))
     call solve_tridiag(a_tri(ks:ke),b_tri(ks:ke),c_tri(ks:ke),d_tri(ks:ke),E_iw(i,j,ks:ke,taup1),ke-ks+1)
    endif
  enddo
 enddo

 !---------------------------------------------------------------------------------
 ! store IW dissipation 
 !---------------------------------------------------------------------------------
 iw_diss = alpha_c*maxE_iw(:,:,:)*E_iw(:,:,:,taup1)

 if (enable_idemix_hor_diffusion) then
 !---------------------------------------------------------------------------------
 ! add tendency due to lateral diffusion
 !---------------------------------------------------------------------------------

  do k=1,nz
   do j=js_pe,je_pe
    do i=is_pe-1,ie_pe
      flux_east(i,j,k)=tau_h*0.5*(v0(i+1,j,k)+v0(i,j,k)) * &
          (v0(i+1,j,k)*E_iw(i+1,j,k,tau)-v0(i,j,k)*E_iw(i,j,k,tau))/(cost(j)*dxu(i))*maskUtr(i,j,k)
    enddo
   enddo
   do j=js_pe-1,je_pe
    do i=is_pe,ie_pe
     flux_north(i,j,k)= tau_h*0.5*(v0(i,j+1,k)+v0(i,j,k)) * &
          (v0(i,j+1,k)*E_iw(i,j+1,k,tau)-v0(i,j,k)*E_iw(i,j,k,tau))/dyu(j)*maskVtr(i,j,k)*cosu(j)
    enddo
   enddo
  enddo

  do j=js_pe,je_pe
    do i=is_pe,ie_pe
     E_iw(i,j,:,taup1)= E_iw(i,j,:,taup1) + dt_tracer*maskW(i,j,:)* &
                                  (( flux_east(i,j,:) - flux_east(i-1,j,:))/(cost(j)*dxt(i)) &
                                  +(flux_north(i,j,:) -flux_north(i,j-1,:))/(cost(j)*dyt(j)) )
   enddo
  enddo
 endif

 if (enable_idemix_hor_diffusion_iter) then
 !---------------------------------------------------------------------------------
 ! add tendency due to lateral diffusion with iterative method in case of high resolution
 !---------------------------------------------------------------------------------

  fxa = dt_tracer/idemix_hor_diffusion_iter
  do n=1,idemix_hor_diffusion_iter

    call border_exchg_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,E_iw(:,:,:,taup1)) 
    call setcyclic_xyz   (is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,E_iw(:,:,:,taup1))
    call set_obc_boundary_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,E_iw(:,:,:,taup1)) 
    do k=1,nz
     do j=js_pe,je_pe
      do i=is_pe-1,ie_pe
        flux_east(i,j,k)=tau_h*0.5*(v0(i+1,j,k)+v0(i,j,k)) * &
            (v0(i+1,j,k)*E_iw(i+1,j,k,taup1)-v0(i,j,k)*E_iw(i,j,k,taup1))/(cost(j)*dxu(i))*maskUtr(i,j,k)
      enddo
     enddo
     do j=js_pe-1,je_pe
      do i=is_pe,ie_pe
       flux_north(i,j,k)= tau_h*0.5*(v0(i,j+1,k)+v0(i,j,k)) * &
            (v0(i,j+1,k)*E_iw(i,j+1,k,taup1)-v0(i,j,k)*E_iw(i,j,k,taup1))/dyu(j)*maskVtr(i,j,k)*cosu(j)
      enddo
     enddo
    enddo

    do j=js_pe,je_pe
      do i=is_pe,ie_pe
       E_iw(i,j,:,taup1)= E_iw(i,j,:,taup1) + fxa*maskW(i,j,:)* &
                                    (( flux_east(i,j,:) - flux_east(i-1,j,:))/(cost(j)*dxt(i)) &
                                    +(flux_north(i,j,:) -flux_north(i,j-1,:))/(cost(j)*dyt(j)) )
     enddo
    enddo
  enddo

 endif

 !---------------------------------------------------------------------------------
 ! add tendency due to advection
 !---------------------------------------------------------------------------------
 if (enable_idemix_superbee_advection) then
  call adv_flux_superbee_wgrid(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,flux_east,flux_north,flux_top,E_iw(:,:,:,tau))
 endif
 if (enable_idemix_dst3_advection) then
  call adv_flux_dst3_wgrid(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,flux_east,flux_north,flux_top,E_iw(:,:,:,tau))
 endif
 if (enable_idemix_superbee_advection .or. enable_idemix_dst3_advection) then
  do j=js_pe,je_pe
   do i=is_pe,ie_pe
      dE_iw(i,j,:,tau)=maskW(i,j,:)* (-( flux_east(i,j,:)-  flux_east(i-1,j,:))/(cost(j)*dxt(i)) &
                                      -(flux_north(i,j,:)- flux_north(i,j-1,:))/(cost(j)*dyt(j)) )
   enddo
  enddo
  k=1; dE_iw(:,:,k,tau)=dE_iw(:,:,k,tau)-flux_top(:,:,k)/dzw(k)   
  do k=2,nz-1
   dE_iw(:,:,k,tau)=dE_iw(:,:,k,tau)-(flux_top(:,:,k)- flux_top(:,:,k-1))/dzw(k)   
  enddo
  k=nz
  dE_iw(:,:,k,tau)=dE_iw(:,:,k,tau)-(flux_top(:,:,k)- flux_top(:,:,k-1))/(0.5*dzw(k))   
  if (enable_idemix_AB_time_stepping ) then 
    !---------------------------------------------------------------------------------
    ! Adam Bashforth time stepping
    !---------------------------------------------------------------------------------
     E_iw(:,:,:,taup1)=E_iw(:,:,:,taup1)+dt_tracer*( (1.5+AB_eps)*dE_iw(:,:,:,tau) &
                          - ( 0.5+AB_eps)*dE_iw(:,:,:,taum1))
  else
     E_iw(:,:,:,taup1)=E_iw(:,:,:,taup1)+dt_tracer*dE_iw(:,:,:,tau)
  endif
  !---------------------------------------------------------------------------------
  ! open boundaries
  !---------------------------------------------------------------------------------
  call set_obc_tracer_wgrid(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,E_iw,.true.) 
 endif

end subroutine integrate_idemix1





function gofx2(x)
!=======================================================================
! a function g(x) 
!=======================================================================
 implicit none
 real*8 :: gofx2,x,c
 real*8, parameter :: pi = 3.14159265358979323846264338327950588
 x=max(3d0,x)
 c= 1.-(2./pi)*asin(1./x)
 gofx2 = 2/pi/c*0.9*x**(-2./3.)*(1-exp(-x/4.3))
end function gofx2

function hofx1(x)
!=======================================================================
! a function h(x) 
!=======================================================================
 implicit none
 real*8 :: hofx1,x
 real*8, parameter :: pi = 3.14159265358979323846264338327950588
 hofx1 = (2./pi)/(1.-(2./pi)*asin(1./x)) * (x-1.)/(x+1.)
end function hofx1

