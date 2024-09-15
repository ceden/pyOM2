




subroutine rossmix2_eddy_advect
 use main_module   
 use rossmix2_module   
 implicit none
 integer :: i,j,k,ks
 real*8 :: dtr(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz)
 real*8 :: aloc(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz),fxa
 
 call rossmix2_eddy_advect_tr_superbee(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,temp(:,:,:,tau),dtr)
 do j=js_pe,je_pe
   do i=is_pe,ie_pe
     temp(i,j,:,taup1)=temp(i,j,:,taup1)+dt_tracer*dtr(i,j,:)
   enddo
 enddo

 if (enable_conserve_energy  .and. .false. ) then
 
  do k=1,nz
   do j=js_pe-onx+1,je_pe+onx-1
    do i=is_pe-onx+1,ie_pe+onx-1
     fxa = int_drhodT(i,j,k,tau)
     aloc(i,j,k) =+0.5*grav/rho_0*( (int_drhodT(i+1,j,k,tau)-fxa)*flux_east(i  ,j,k) &
                                   +(fxa-int_drhodT(i-1,j,k,tau))*flux_east(i-1,j,k) ) /(dxt(i)*cost(j))  &
                  +0.5*grav/rho_0*( (int_drhodT(i,j+1,k,tau)-fxa)*flux_north(i,j  ,k) &
                                   +(fxa-int_drhodT(i,j-1,k,tau))*flux_north(i,j-1,k) ) /(dyt(j)*cost(j)) 
    enddo
   enddo
  end do

  do j=js_pe-onx,je_pe+onx
   do i=is_pe-onx,ie_pe+onx
     ks=kbot(i,j)
     if (ks>0) then
      k=ks; K_diss_gm(i,j,k) =  &
                         0.5*(aloc(i,j,k)+aloc(i,j,k+1)) + 0.5*aloc(i,j,k)*dzw(max(1,k-1))/dzw(k)
      do k=ks+1,nz-1
       K_diss_gm(i,j,k) = 0.5*(aloc(i,j,k) +aloc(i,j,k+1))
      enddo
      k=nz; K_diss_gm(i,j,k) = aloc(i,j,k)
     endif
   enddo
  enddo

  do k=1,nz-1
   do j=js_pe,je_pe
    do i=is_pe,ie_pe
     fxa = (-int_drhodT(i,j,k+1,tau) +int_drhodT(i,j,k,tau))/dzw(k)
     K_diss_gm(i,j,k) = K_diss_gm(i,j,k)  -grav/rho_0*fxa*flux_top(i,j,k)*maskW(i,j,k)  
    enddo
   enddo
  end do
 endif



 call rossmix2_eddy_advect_tr_superbee(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,salt(:,:,:,tau),dtr)
 do j=js_pe,je_pe
   do i=is_pe,ie_pe
     salt(i,j,:,taup1)=salt(i,j,:,taup1)+dt_tracer*dtr(i,j,:)
   enddo
 enddo

 if (enable_conserve_energy .and. .false. ) then

  do k=1,nz
   do j=js_pe-onx+1,je_pe+onx-1
    do i=is_pe-onx+1,ie_pe+onx-1
     fxa = int_drhodS(i,j,k,tau)
     aloc(i,j,k) =+0.5*grav/rho_0*( (int_drhodS(i+1,j,k,tau)-fxa)*flux_east(i  ,j,k) &
                                   +(fxa-int_drhodS(i-1,j,k,tau))*flux_east(i-1,j,k) ) /(dxt(i)*cost(j))  &
                  +0.5*grav/rho_0*( (int_drhodS(i,j+1,k,tau)-fxa)*flux_north(i,j  ,k) &
                                   +(fxa-int_drhodS(i,j-1,k,tau))*flux_north(i,j-1,k) ) /(dyt(j)*cost(j)) 
    enddo
   enddo
  end do

  do j=js_pe-onx,je_pe+onx
   do i=is_pe-onx,ie_pe+onx
     ks=kbot(i,j)
     if (ks>0) then
      k=ks; K_diss_gm(i,j,k) = K_diss_gm(i,j,k)+ &
                         0.5*(aloc(i,j,k)+aloc(i,j,k+1)) + 0.5*aloc(i,j,k)*dzw(max(1,k-1))/dzw(k)
      do k=ks+1,nz-1
       K_diss_gm(i,j,k) = K_diss_gm(i,j,k)+ 0.5*(aloc(i,j,k) +aloc(i,j,k+1))
      enddo
      k=nz; K_diss_gm(i,j,k) = K_diss_gm(i,j,k)+ aloc(i,j,k)
     endif
   enddo
  enddo

  do k=1,nz-1
   do j=js_pe,je_pe
    do i=is_pe,ie_pe
     fxa = (-int_drhodS(i,j,k+1,tau) +int_drhodS(i,j,k,tau))/dzw(k)
     K_diss_gm(i,j,k) = K_diss_gm(i,j,k)  -grav/rho_0*fxa*flux_top(i,j,k)*maskW(i,j,k)  
    enddo
   enddo
  end do

 endif
end subroutine rossmix2_eddy_advect




subroutine rossmix2_eddy_advect_tr_2nd(is_,ie_,js_,je_,nz_,var,dvar)
 use main_module   
 use rossmix2_module   
 implicit none
 integer :: is_,ie_,js_,je_,nz_
 real*8 :: var(is_:ie_,js_:je_,nz_),dvar(is_:ie_,js_:je_,nz_)
 integer :: i,j,k

 do i=is_pe-1,ie_pe
   flux_east(i,:,:) = 0.5*(var(i,:,:)+var(i+1,:,:))*ue(i,:,:)*maskU(i,:,:)
 enddo
 do j=js_pe-1,je_pe
   flux_north(:,j,:) = cosu(j)*0.5*(var(:,j,:)+var(:,j+1,:))*ve(:,j,:)*maskV(:,j,:)
 enddo
 do k=1,nz-1
   flux_top(:,:,k)=0.5*(var(:,:,k)+var(:,:,k+1))*we(:,:,k)*maskW(:,:,k)
 enddo
 flux_top(:,:,nz)=0.0
 
 do j=js_pe,je_pe
   do i=is_pe,ie_pe
      dvar(i,j,:)=maskT(i,j,:)* (-( flux_east(i,j,:)-  flux_east(i-1,j,:))/(cost(j)*dxt(i)) &
                                  -(flux_north(i,j,:)- flux_north(i,j-1,:))/(cost(j)*dyt(j)) )
   enddo
 enddo
 k=1; dvar(:,:,k)=dvar(:,:,k)-maskT(:,:,k)*flux_top(:,:,k)/dzt(k)   
 do k=2,nz
   dvar(:,:,k)=dvar(:,:,k)-maskT(:,:,k)*(flux_top(:,:,k)- flux_top(:,:,k-1))/dzt(k)   
 enddo
end subroutine rossmix2_eddy_advect_tr_2nd




subroutine rossmix2_eddy_advect_tr_superbee(is_,ie_,js_,je_,nz_,var,dvar)
 use main_module   
 use rossmix2_module   
 implicit none
 integer :: is_,ie_,js_,je_,nz_
 real*8 :: var(is_:ie_,js_:je_,nz_),dvar(is_:ie_,js_:je_,nz_)
 integer :: i,j,k,km1,kp2
 real*8 :: Rjp,Rj,Rjm,uCFL,Cr,Limiter
 Limiter(Cr)=max(0.D0,max(min(1.D0,2.D0*Cr), min(2.D0,Cr))) 

 do k=1,nz
   do j=js_pe,je_pe
        do i=is_pe-1,ie_pe
         uCFL = ABS( ue(i,j,k)*dt_tracer/(cost(j)*dxt(i)) )
         Rjp=(var(i+2,j,k)-var(i+1,j,k))*maskU(i+1,j,k)
         Rj =(var(i+1,j,k)-var(i  ,j,k))*maskU(i  ,j,k)
         Rjm=(var(i  ,j,k)-var(i-1,j,k))*maskU(i-1,j,k)
         IF (Rj.NE.0.) THEN
          IF (ue(i,j,k).GT.0) THEN; Cr=Rjm/Rj; ELSE; Cr=Rjp/Rj; ENDIF
         ELSE
          IF (ue(i,j,k).GT.0) THEN; Cr=Rjm*1.E20; ELSE; Cr=Rjp*1.E20; ENDIF
         ENDIF
         Cr=Limiter(Cr)
         flux_east(i,j,k) = (ue(i,j,k)*(var(i+1,j,k)+var(i,j,k))*0.5d0   &
                                -ABS(ue(i,j,k))*((1.-Cr)+uCFL*Cr)*Rj*0.5d0)*maskU(i,j,k)
        enddo
   enddo
 enddo

 do k=1,nz
   do j=js_pe-1,je_pe
        do i=is_pe,ie_pe
         Rjp=(var(i,j+2,k)-var(i,j+1,k))*maskV(i,j+1,k)
         Rj =(var(i,j+1,k)-var(i,j  ,k))*maskV(i,j  ,k)
         Rjm=(var(i,j  ,k)-var(i,j-1,k))*maskV(i,j-1,k)
         uCFL = ABS( cosu(j)*ve(i,j,k)*dt_tracer/(cost(j)*dyt(j)) )
         IF (Rj.NE.0.) THEN
          IF (ve(i,j,k).GT.0) THEN; Cr=Rjm/Rj; ELSE; Cr=Rjp/Rj; ENDIF
         ELSE
          IF (ve(i,j,k).GT.0) THEN; Cr=Rjm*1.E20; ELSE; Cr=Rjp*1.E20; ENDIF
         ENDIF
         Cr=Limiter(Cr)
         flux_north(i,j,k) = (cosu(j)*ve(i,j,k)*(var(i,j+1,k)+var(i,j,k))*0.5d0   &
                    -ABS(cosu(j)*ve(i,j,k))*((1.-Cr)+uCFL*Cr)*Rj*0.5d0)*maskV(i,j,k)
        enddo
   enddo
 enddo
 
 do k=1,nz-1
   kp2=min(nz,k+2); !if (kp2>nphi) kp2=3
   km1=max(1,k-1) !if (km1<1) km1=nphi-2
   do j=js_pe,je_pe
      do i=is_pe,ie_pe
         Rjp=(var(i,j,kp2)-var(i,j,k+1))*maskW(i,j,k+1)
         Rj =(var(i,j,k+1)-var(i,j,k  ))*maskW(i,j,k  )
         Rjm=(var(i,j,k  )-var(i,j,km1))*maskW(i,j,km1)
         uCFL = ABS( we(i,j,k)*dt_tracer/dzt(k) )
         IF (Rj.NE.0.) THEN
          IF (we(i,j,k).GT.0) THEN; Cr=Rjm/Rj; ELSE; Cr=Rjp/Rj; ENDIF
         ELSE
          IF (we(i,j,k).GT.0) THEN; Cr=Rjm*1.E20; ELSE; Cr=Rjp*1.E20; ENDIF
         ENDIF
         Cr=Limiter(Cr)
         flux_top(i,j,k) = (we(i,j,k)*(var(i,j,k+1)+var(i,j,k))*0.5d0   &
                                -ABS(we(i,j,k))*((1.-Cr)+uCFL*Cr)*Rj*0.5d0)*maskW(i,j,k)
      enddo
   enddo
 enddo
 flux_top(:,:,nz)=0.0
 
 do j=js_pe,je_pe
   do i=is_pe,ie_pe
      dvar(i,j,:)=maskT(i,j,:)* (-( flux_east(i,j,:)-  flux_east(i-1,j,:))/(cost(j)*dxt(i)) &
                                  -(flux_north(i,j,:)- flux_north(i,j-1,:))/(cost(j)*dyt(j)) )
   enddo
 enddo
 k=1; dvar(:,:,k)=dvar(:,:,k)-maskT(:,:,k)*flux_top(:,:,k)/dzt(k)   
 do k=2,nz
   dvar(:,:,k)=dvar(:,:,k)-maskT(:,:,k)*(flux_top(:,:,k)- flux_top(:,:,k-1))/dzt(k)   
 enddo
end subroutine rossmix2_eddy_advect_tr_superbee






subroutine rossmix2_friction 
use main_module   
 use rossmix2_module   
 implicit none
 integer :: i,j,k,n
 real*8 :: stress_uz(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz)
 real*8 :: stress_vz(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz)
 !real*8 :: diss(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz)
 real*8 :: aloc(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx)
 real*8 :: ediss(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx)
 real*8 :: uloc(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx)
 real*8 :: vloc(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx)
 real*8 :: r_fac(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx)
 
if (.false.) then 
 ! ubar = int U dz/h
 uloc=0; aloc=0
 do k=1,nz
   uloc = uloc + u(:,:,k,tau)*maskU(:,:,k)*dzt(k)
   aloc = aloc + maskU(:,:,k)*dzt(k)
 enddo
 where (aloc/=0) uloc=uloc/aloc
 
 ! vbar = int v dz/h
 vloc=0; aloc=0
 do k=1,nz
   vloc = vloc + v(:,:,k,tau)*maskV(:,:,k)*dzt(k)
   aloc = aloc + maskV(:,:,k)*dzt(k)
 enddo
 where (aloc/=0) vloc=vloc/aloc
 
 ! int gamma e^2 dalpha
 ediss=0
 do n=2,nphi-1
  ediss = ediss + E_r(:,:,n,taup1)*E_r(:,:,n,tau)*rossmix2_damping*maskTp(:,:,n)*dphit(n)
 enddo
 
 ! \int (u u_bar + v v_bar) dz
 aloc = 0
 do k=1,nz
  aloc = aloc + (u(:,:,k,tau)*uloc(:,:)*maskU(:,:,k) + v(:,:,k,tau)*vloc(:,:)*maskV(:,:,k))*dzt(k)
 enddo
 ! r = 0.5 int gamma e^2 dalpha / \int (u u_bar + v v_bar) dz
 aloc = sign( max(1e-12,abs(aloc)), aloc )
 r_fac=0.75*ediss/aloc 
 
 ! + r u_bar
 do k=1,nz
  du_mix(:,:,k)= du_mix(:,:,k) + maskU(:,:,k)*r_fac*uloc 
 enddo
 
 do k=1,nz
    dv_mix(:,:,k)= dv_mix(:,:,k) + maskV(:,:,k)*r_fac*vloc
 enddo
endif 
 
 !K_diss_ross = 0.
 !---------------------------------------------------------------------------------
 ! divergence of horizontal eddy zonal momentum flux 
 !---------------------------------------------------------------------------------
 flux_east(:,:,:) = -stress_ux(:,:,:)
 do j=js_pe-1,je_pe
   flux_north(:,j,:)= -stress_uy(:,j,:)*cosu(j)
 enddo
 do j=js_pe,je_pe
   do i=is_pe,ie_pe
    du_mix(i,j,:)= du_mix(i,j,:) + maskU(i,j,:)*((flux_east(i,j,:) -   flux_east(i-1,j,:))/(cost(j)*dxu(i)) &
                                                +(flux_north(i,j,:) - flux_north(i,j-1,:))/(cost(j)*dyt(j)) )
  enddo
 enddo
 !---------------------------------------------------------------------------------
 ! diagnose dissipation by that flux
 !---------------------------------------------------------------------------------
  !do k=1,nz
  ! do j=js_pe,je_pe
  !  do i=is_pe-1,ie_pe
  !   diss(i,j,k) =0.5*((u(i+1,j,k,tau)-u(i,j,k,tau))*flux_east(i,j,k) &
  !                    +(u(i,j,k,tau)-u(i-1,j,k,tau))*flux_east(i-1,j,k))/(cost(j)*dxu(i))  &
  !               +0.5*((u(i,j+1,k,tau)-u(i,j,k,tau))*flux_north(i,j,k)+ &
  !                     (u(i,j,k,tau)-u(i,j-1,k,tau))*flux_north(i,j-1,k))/(cost(j)*dyt(j)) 
  !  enddo
  ! enddo
  !enddo
  !call calc_diss(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,diss,K_diss_ross,'U')

 !---------------------------------------------------------------------------------
 ! divergence of horizontal eddy meridional momentum flux 
 !---------------------------------------------------------------------------------
 flux_east(:,:,:) = -stress_vx(:,:,:)
 do j=js_pe-1,je_pe
   flux_north(:,j,:)= -stress_vy(:,j,:)*cost(j+1)
 enddo
 do j=js_pe,je_pe
   do i=is_pe,ie_pe
    dv_mix(i,j,:)= dv_mix(i,j,:) + maskV(i,j,:)*( (flux_east(i,j,:) -   flux_east(i-1,j,:))/(cosu(j)*dxt(i))  &
                                                 +(flux_north(i,j,:) - flux_north(i,j-1,:))/(dyu(j)*cosu(j)) )
  enddo
 enddo

 !---------------------------------------------------------------------------------
 ! diagnose dissipation by that flux
 !---------------------------------------------------------------------------------
 ! do k=1,nz
 !  do j=js_pe-1,je_pe
 !   do i=is_pe,ie_pe
 !    diss(i,j,k) =0.5*((v(i+1,j,k,tau)-v(i,j,k,tau))*flux_east(i,j,k)+ &
 !                      (v(i,j,k,tau)-v(i-1,j,k,tau))*flux_east(i-1,j,k))/(cosu(j)*dxt(i)) &
 !               + 0.5*((v(i,j+1,k,tau)-v(i,j,k,tau))*flux_north(i,j,k)+ &
 !                      (v(i,j,k,tau)-v(i,j-1,k,tau))*flux_north(i,j-1,k))/(cosu(j)*dyu(j)) 
 !   enddo
 !  enddo
 ! enddo
  !call calc_diss(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,diss,K_diss_ross,'V')



 return

 
 do k=1,nz-1
    stress_uz(:,:,k) = -coriolis_t(:,:)*Bx(:,:,k)
    stress_vz(:,:,k) = -coriolis_t(:,:)*By(:,:,k)
 enddo
 call tgrid_to_ugrid( is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,stress_uz,stress_uz)  
 call tgrid_to_vgrid( is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,stress_vz,stress_vz)  

 ! no flux from bottom or surface
 do k=1,nz-1
   stress_uz(:,:,k) = stress_uz(:,:,k)*maskU(:,:,k)*maskU(:,:,k+1)
   stress_vz(:,:,k) = stress_vz(:,:,k)*maskV(:,:,k)*maskV(:,:,k+1)
 enddo
 stress_uz(:,:,nz)=0; stress_vz(:,:,nz)=0;

 !---------------------------------------------------------------------------------
 ! divergence of vertical eddy zonal momentum flux 
 !---------------------------------------------------------------------------------
 do k=1,nz-1
  do j=js_pe,je_pe
   do i=is_pe,ie_pe
    flux_top(i,j,k)=  stress_uz(i,j,k)
   enddo
  enddo 
 enddo
 flux_top(:,:,nz) = 0.0 
 k=1; du_mix(:,:,k) = du_mix(:,:,k) + flux_top(:,:,k)/dzt(k)*maskU(:,:,k)
 do k=2,nz
   du_mix(:,:,k) = du_mix(:,:,k) + (flux_top(:,:,k)-flux_top(:,:,k-1))/dzt(k)*maskU(:,:,k)
 enddo

 !---------------------------------------------------------------------------------
 ! divergence of vertical eddy meridional momentum flux 
 !---------------------------------------------------------------------------------
 do k=1,nz-1
  do j=js_pe,je_pe
   do i=is_pe,ie_pe
    flux_top(i,j,k)=  stress_vz(i,j,k)
   enddo
  enddo     
 enddo
 flux_top(:,:,nz) = 0.0 
 k=1; dv_mix(:,:,k) = dv_mix(:,:,k) + flux_top(:,:,k)/dzt(k)*maskV(:,:,k)
 do k=2,nz
   dv_mix(:,:,k) = dv_mix(:,:,k) + (flux_top(:,:,k)-flux_top(:,:,k-1))/dzt(k)*maskV(:,:,k)
 enddo
 
end subroutine rossmix2_friction
