




subroutine idemix_get_delta(ltau)
!=======================================================================
!
!=======================================================================
 use main_module   
 use idemix_module   
 implicit none
 integer :: ltau
 integer :: k,i,j,ks

 !---------------------------------------------------------------------------------
 ! delta F = - tau (c0 F)_z  - tau Tc F
 !---------------------------------------------------------------------------------
 do j=js_pe-onx,je_pe+onx
  do i=is_pe-onx,ie_pe+onx-1
    ks=max(kbot(i,j),kbot(i+1,j))
    if (ks>0) then
      do k=ks,nz-1
        flux_top(i,j,k) = 0.5*(c0_u(i,j,k)+c0_u(i,j,k+1))*F_s(i,j,k+1,ltau )*maskU(i,j,k+1)
      enddo
      k=ks;F_d(i,j,k,ltau ) = - tau_v*( flux_top(i,j,k) - c0_u(i,j,k)*F_s(i,j,k,ltau ) )/dzw(k) &
                              - tau_v*Tc_u(i,j,k)*0.5*(F_s(i,j,k,ltau )+ F_s(i,j,k+1,ltau ))
      do k=ks+1,nz-1
        F_d(i,j,k,ltau )   = - tau_v*( flux_top(i,j,k)-flux_top(i,j,k-1) )/dzw(k)  &
                             - tau_v*Tc_u(i,j,k)*0.5*(F_s(i,j,k,ltau )+ F_s(i,j,k+1,ltau ))
      enddo
      !k=nz; F_d(i,j,k,ltau ) = 0.25*forc_iw_surface_u(i,j)/c0_u(i,j,k-1)
    endif
  enddo
 enddo

 !---------------------------------------------------------------------------------
 ! delta G = - tau (c0 G)_z  + tau Tc G
 !---------------------------------------------------------------------------------
 do j=js_pe-onx,je_pe+onx
  do i=is_pe-onx,ie_pe+onx-1
    ks=max(kbot(i,j),kbot(i+1,j))
    if (ks>0) then
      do k=ks,nz-1
        flux_top(i,j,k) = 0.5*(c0_u(i,j,k)+c0_u(i,j,k+1))*G_s(i,j,k+1,ltau )*maskU(i,j,k+1)
      enddo
      k=ks;G_d(i,j,k,ltau ) = - tau_v*( flux_top(i,j,k) - c0_u(i,j,k)*G_s(i,j,k,ltau ) )/dzw(k) &
                              + tau_v*Tc_u(i,j,k)*0.5*(G_s(i,j,k,ltau )+ G_s(i,j,k+1,ltau ))
      do k=ks+1,nz-1
        G_d(i,j,k,ltau )   = - tau_v*( flux_top(i,j,k)-flux_top(i,j,k-1) )/dzw(k)  &
                             + tau_v*Tc_u(i,j,k)*0.5*(G_s(i,j,k,ltau )+ G_s(i,j,k+1,ltau ))
      enddo
      !k=nz; G_d(i,j,k,ltau ) = 0.25*forc_iw_surface_u(i,j)/c0_u(i,j,k-1)
    endif
  enddo
 enddo

 !---------------------------------------------------------------------------------
 ! delta O = - tau (c0 O)_z  - tau Tc O
 !---------------------------------------------------------------------------------
 do j=js_pe-onx,je_pe+onx-1
  do i=is_pe-onx,ie_pe+onx
    ks=max(kbot(i,j),kbot(i,j+1))
    if (ks>0) then
      do k=ks,nz-1
        flux_top(i,j,k) = 0.5*(c0_v(i,j,k)+c0_v(i,j,k+1))*O_s(i,j,k+1,ltau )*maskV(i,j,k+1)
      enddo
      k=ks;O_d(i,j,k,ltau ) = - tau_v*( flux_top(i,j,k) - c0_v(i,j,k)*O_s(i,j,k,ltau ) )/dzw(k) &
                              - tau_v*Tc_v(i,j,k)*0.5*(O_s(i,j,k,ltau )+ O_s(i,j,k+1,ltau ))
      do k=ks+1,nz-1
        O_d(i,j,k,ltau )   = - tau_v*( flux_top(i,j,k)-flux_top(i,j,k-1) )/dzw(k)  &
                             - tau_v*Tc_v(i,j,k)*0.5*(O_s(i,j,k,ltau )+ O_s(i,j,k+1,ltau ))
      enddo
      !k=nz; O_d(i,j,k,ltau ) = 0.25*forc_iw_surface_v(i,j)/c0_v(i,j,k-1)
    endif
  enddo
 enddo

 !---------------------------------------------------------------------------------
 ! delta P = - tau (c0 P)_z  + tau Tc P
 !---------------------------------------------------------------------------------
 do j=js_pe-onx,je_pe+onx-1
  do i=is_pe-onx,ie_pe+onx
    ks=max(kbot(i,j),kbot(i,j+1))
    if (ks>0) then
      do k=ks,nz-1
        flux_top(i,j,k) = 0.5*(c0_v(i,j,k)+c0_v(i,j,k+1))*P_s(i,j,k+1,ltau )*maskV(i,j,k+1)
      enddo
      k=ks;P_d(i,j,k,ltau ) = - tau_v*( flux_top(i,j,k) - c0_v(i,j,k)*P_s(i,j,k,ltau ) )/dzw(k) &
                              + tau_v*Tc_v(i,j,k)*0.5*(P_s(i,j,k,ltau )+ P_s(i,j,k+1,ltau ))
      do k=ks+1,nz-1
        P_d(i,j,k,ltau )   = - tau_v*( flux_top(i,j,k)-flux_top(i,j,k-1) )/dzw(k)  &
                             + tau_v*Tc_v(i,j,k)*0.5*(P_s(i,j,k,ltau )+ P_s(i,j,k+1,ltau ))
      enddo
      !k=nz; P_d(i,j,k,ltau ) = 0.25*forc_iw_surface_v(i,j)/c0_v(i,j,k-1)
    endif
  enddo
 enddo

end subroutine idemix_get_delta








subroutine idemix_friction(ltau)
!=======================================================================
!  wave-induced gravity wave drag from idemix
!=======================================================================
 use main_module   
 use idemix_module   
 implicit none
 integer :: i,j,k
 integer :: ltau
 real*8 :: diss(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz)

 !---------------------------------------------------------------------------------
 ! explicit vertical friction of zonal momentum
 ! F_u = sqrt(8)/pi Omega (Delta E_east - Delta E_west)
 !     = sqrt(8)/pi Omega (F_d - G_d)
 ! u_t = ...  - (F_u)_z
 !---------------------------------------------------------------------------------
 do k=1,nz-1
   do i=is_pe-1,ie_pe
    flux_top(i,:,k)=  sqrt(8.)/pi*Om_id3_u(i,:,k)* &
                   (G_d(i,:,k,ltau)-F_d(i,:,k,ltau))*maskU(i,:,k+1)*maskU(i,:,k)
   enddo
 enddo
 flux_top(:,:,nz) = 0.0 
 k=1; du_mix(:,:,k) = du_mix(:,:,k)+flux_top(:,:,k)/dzt(k)*maskU(:,:,k)
 do k=2,nz
   du_mix(:,:,k) = du_mix(:,:,k)+(flux_top(:,:,k)-flux_top(:,:,k-1))/dzt(k)*maskU(:,:,k)
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
 call wgrid_to_tgrid(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,diss,diss)
 call ugrid_to_tgrid(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,diss,diss)
 K_diss_idemix = diss

 !---------------------------------------------------------------------------------
 ! explicit vertical friction of meridional momentum
 ! F_v = sqrt(8)/pi Omega (Delta E_north - Delta E_south)
 !     = sqrt(8)/pi Omega (O_d - P_d)
 ! v_t = ...  - (F_v)_z
 !---------------------------------------------------------------------------------
 do k=1,nz-1
  do j=js_pe-1,je_pe
    flux_top(:,j,k)=  sqrt(8.)/pi*Om_id3_v(:,j,k)* &
                   (P_d(:,j,k,ltau)-O_d(:,j,k,ltau))*maskV(:,j,k+1)*maskV(:,j,k)
   enddo
 enddo
 flux_top(:,:,nz) = 0.0 
 k=1; dv_mix(:,:,k) = dv_mix(:,:,k)+ flux_top(:,:,k)/dzt(k)*maskV(:,:,k)
 do k=2,nz
   dv_mix(:,:,k) = dv_mix(:,:,k)+ (flux_top(:,:,k)-flux_top(:,:,k-1))/dzt(k)*maskV(:,:,k)
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
 call wgrid_to_tgrid(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,diss,diss)
 call vgrid_to_tgrid(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,diss,diss)
 K_diss_idemix = K_diss_idemix + diss

 !---------------------------------------------------------------------------------
 ! Diagnose viscosities
 !---------------------------------------------------------------------------------
 where (Tc_u /= 0d0) 
      viscU_idemix = maskU*(sqrt(8.)/pi*Om_id3_u)**2*(G_d(:,:,:,ltau)-F_d(:,:,:,ltau))/Tc_u 
 elsewhere
      viscU_idemix = 0d0
 endwhere

 where (Tc_v /= 0d0) 
      viscV_idemix = maskV*(sqrt(8.)/pi*Om_id3_v)**2*(P_d(:,:,:,ltau)-O_d(:,:,:,ltau))/Tc_v 
 elsewhere
      viscV_idemix = 0d0
 endwhere

end subroutine idemix_friction






