
subroutine biharmonic_thickness_skew_tensor
 use main_module 
 use biharmonic_thickness_module
 use diagnostics_module
 implicit none  
 integer :: i,j,k,ip,jp,kr,n
 real*8 :: get_drhodT,get_drhodS 
 real*8 :: Ai_ez(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz,0:1,0:1) 
 real*8 :: Ai_nz(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz,0:1,0:1) 
 real*8 :: Ai_bx(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz,0:1,0:1)
 real*8 :: Ai_by(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz,0:1,0:1) 
 real*8 :: drdTS(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz,2)
 real*8 :: ddzt(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz,2)
 real*8 :: ddxt(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz,2)
 real*8 :: ddyt(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz,2)
 real*8 :: ddx2t(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz,2)
 real*8 :: ddy2t(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz,2)
 real*8 :: diff_loc(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz)
 real*8 :: drodx2e,drodze,drody2n,drodzn,drodx2b,drody2b,drodzb
 real*8 :: taper,sxe,syn,sxb,syb,clip,z1,z2
 real*8,parameter :: epsln=1.D-20 
 
 taper(sxe,z1,z2) = 0.5*(1.+tanh((-abs(sxe)+biha_slopec)/biha_dslope)) *0.5*(1.+tanh((-z1-z2)/10.)) 
 !taper(sxe) = 0.5*(1.+tanh((-abs(sxe)+biha_slopec)/biha_dslope)) 
!-----------------------------------------------------------------------
 ! statement functions for density triads
!-----------------------------------------------------------------------
 drodx2e(i,j,k,ip)    = drdTS(i+ip,j,k,1)*ddx2t(i,j,k,1)         + drdTS(i+ip,j,k,2)*ddx2t(i,j,k,2) 
 drodze(i,j,k,ip,kr)  = drdTS(i+ip,j,k,1)*ddzt(i+ip,j,k+kr-1,1)  + drdTS(i+ip,j,k,2)*ddzt(i+ip,j,k+kr-1,2)
 drody2n(i,j,k,jp)    = drdTS(i,j+jp,k,1)*ddy2t(i,j,k,1)         + drdTS(i,j+jp,k,2)*ddy2t(i,j,k,2) 
 drodzn(i,j,k,jp,kr)  = drdTS(i,j+jp,k,1)*ddzt(i,j+jp,k+kr-1,1)  + drdTS(i,j+jp,k,2)*ddzt(i,j+jp,k+kr-1,2)
 drodx2b(i,j,k,ip,kr) = drdTS(i,j,k+kr,1)*ddx2t(i-1+ip,j,k+kr,1) + drdTS(i,j,k+kr,2)*ddx2t(i-1+ip,j,k+kr,2) 
 drody2b(i,j,k,jp,kr) = drdTS(i,j,k+kr,1)*ddy2t(i,j-1+jp,k+kr,1) + drdTS(i,j,k+kr,2)*ddy2t(i,j-1+jp,k+kr,2) 
 drodzb(i,j,k,kr)     = drdTS(i,j,k+kr,1)*ddzt(i,j,k,1)          + drdTS(i,j,k+kr,2)*ddzt(i,j,k,2)
  
!-----------------------------------------------------------------------
!     drho_dt and drho_ds at centers of T cells
!-----------------------------------------------------------------------
 do k=1,nz
  do j=js_pe-onx,je_pe+onx
   do i=is_pe-onx,ie_pe+onx
    drdTS(i,j,k,1) = get_drhodT(salt(i,j,k,tau),temp(i,j,k,tau),abs(zt(k)))*maskT(i,j,k)
    drdTS(i,j,k,2) = get_drhodS(salt(i,j,k,tau),temp(i,j,k,tau),abs(zt(k)))*maskT(i,j,k)
   enddo
  enddo
 enddo
!-----------------------------------------------------------------------
!     gradients at top face of T cells
!-----------------------------------------------------------------------
 do k=1,nz-1
  do j=js_pe-onx,je_pe+onx
   do i=is_pe-onx,ie_pe+onx
    ddzt(i,j,k,1) = maskW(i,j,k)* (temp(i,j,k+1,tau) - temp(i,j,k,tau))/dzw(k)  
    ddzt(i,j,k,2) = maskW(i,j,k)* (salt(i,j,k+1,tau) - salt(i,j,k,tau))/dzw(k) 
   enddo
  enddo
 enddo
 ddzt(:,:,nz,:)=0
!-----------------------------------------------------------------------
!    T/S gradients at eastern face of T cells
!-----------------------------------------------------------------------
 do j=js_pe-onx,je_pe+onx
  do i=is_pe-onx,ie_pe+onx-1
     ddxt(i,j,:,1) = maskU(i,j,:)*(temp(i+1,j,:,tau)-temp(i,j,:,tau))/(dxu(i)*cost(j))
     ddxt(i,j,:,2) = maskU(i,j,:)*(salt(i+1,j,:,tau)-salt(i,j,:,tau))/(dxu(i)*cost(j))
  enddo
 enddo
 ddxt(ie_pe+onx,:,:,:)=0d0

!-----------------------------------------------------------------------
!     laplace of T/S gradients at eastern face of T cells
!-----------------------------------------------------------------------
 do n=1,2
  do j=js_pe-1,je_pe+1   
   do i=is_pe-2,ie_pe+1
    flux_east(i,j,:) = (ddxt(i+1,j,:,n)-ddxt(i,j,:,n)) /(cost(j)*dxt(i+1))*maskU(i+1,j,:)*maskU(i,j,:)
   enddo
  enddo
  do j=js_pe-2,je_pe+1
   do i=is_pe-1,ie_pe+1
    flux_north(i,j,:) = (ddxt(i,j+1,:,n)-ddxt(i,j,:,n))/dyu(j)*maskU(i,j+1,:)*maskU(i,j,:)*cosu(j)
   enddo 
  enddo 
!  do j=js_pe-1,je_pe+1
!   do i=is_pe-1,ie_pe+1
  do j=js_pe,je_pe
   do i=is_pe-1,ie_pe   
    ddx2t(i,j,:,n)= maskU(i,j,:)*((flux_east(i,j,:) - flux_east(i-1,j,:))/(cost(j)*dxu(i)) &
                                +(flux_north(i,j,:) - flux_north(i,j-1,:))/(cost(j)*dyt(j)) )
   enddo
  enddo 
 enddo
 
!-----------------------------------------------------------------------
!     gradients at northern face of T cells
!-----------------------------------------------------------------------
 do j=js_pe-onx,je_pe+onx-1
  do i=is_pe-onx,ie_pe+onx
    ddyt(i,j,:,1) = maskV(i,j,:)*(temp(i,j+1,:,tau) - temp(i,j,:,tau)) /dyu(j)
    ddyt(i,j,:,2) = maskV(i,j,:)*(salt(i,j+1,:,tau) - salt(i,j,:,tau)) /dyu(j)
  enddo
 enddo
 ddyt(:,je_pe+onx,:,:)=0d0
 
 do n=1,2
  do j=js_pe-1,je_pe+1
   do i=is_pe-2,ie_pe+1
    flux_east(i,j,:) = (ddyt(i+1,j,:,n)-ddyt(i,j,:,n))/(cosu(j)*dxu(i)) *maskV(i+1,j,:)*maskV(i,j,:)
   enddo
  enddo
  do j=js_pe-2,je_pe+1
   do i=is_pe-1,ie_pe+1
    flux_north(i,j,:) = (ddyt(i,j+1,:,n)-ddyt(i,j,:,n) )/dyt(j+1)*cost(j+1)*maskV(i,j,:)*maskV(i,j+1,:)
   enddo
  enddo 
 ! do j=js_pe-1,je_pe+1
 !  do i=is_pe-1,ie_pe+1
  do j=js_pe-1,je_pe
   do i=is_pe,ie_pe   
    ddy2t(i,j,:,n)=  maskV(i,j,:)*( (flux_east(i,j,:) - flux_east(i-1,j,:))/(cosu(j)*dxt(i))  &
                                  +(flux_north(i,j,:) - flux_north(i,j-1,:))/(dyu(j)*cosu(j)) )
   enddo
  enddo
 enddo 
!-----------------------------------------------------------------------
!     Compute Ai_ez on center of east face of T cell.
!-----------------------------------------------------------------------
 do k=2,nz
  do j=js_pe,je_pe
   do i=is_pe-1,ie_pe
    do kr=0,1
     do ip=0,1
       sxe  = -drodx2e(i,j,k,ip)/(min(0d0,drodze(i,j,k,ip,kr))-epsln)  ! i+1, k-1
       clip = taper(sxe,zw(k),mldu(i,j)+mld0_thk)
       Ai_ez(i,j,k,ip,kr) =  clip*sxe*maskU(i,j,k)
     enddo
    enddo
   enddo
  enddo
 enddo
 k=1
 do j=js_pe,je_pe
   do i=is_pe-1,ie_pe
    kr=1
    do ip=0,1
       sxe  = -drodx2e(i,j,k,ip)/(min(0d0,drodze(i,j,k,ip,kr))-epsln)
       clip = taper(sxe,zw(k),mldu(i,j)+mld0_thk)
       Ai_ez(i,j,k,ip,kr) =  clip*sxe*maskU(i,j,k)
    enddo
   enddo
 enddo
 Ai_ez(:,:,1,:,0) =0d0
!-----------------------------------------------------------------------
!     Compute Ai_nz on center of north face of T cell.
!-----------------------------------------------------------------------
 do k=2,nz
  do j=js_pe-1,je_pe
   do i=is_pe,ie_pe
    do kr=0,1
      do jp=0,1
         syn = -drody2n(i,j,k,jp)/(min(0d0,drodzn(i,j,k,jp,kr))-epsln)
         clip = taper(syn,zw(k),mldv(i,j)+mld0_thk)
         Ai_nz(i,j,k,jp,kr) = clip*syn*maskV(i,j,k)
      enddo
     enddo
    enddo
  enddo
 enddo
 k=1
 do j=js_pe-1,je_pe
   do i=is_pe,ie_pe
     kr=1
     do jp=0,1
         syn = -drody2n(i,j,k,jp)/(min(0d0,drodzn(i,j,k,jp,kr))-epsln)
         clip = taper(syn,zw(k),mldv(i,j)+mld0_thk)
         Ai_nz(i,j,k,jp,kr) = clip*syn*maskV(i,j,k)         
     enddo
   enddo
 enddo
 Ai_nz(:,:,1,:,0) = 0d0 
!-----------------------------------------------------------------------
! compute Ai_bx, Ai_by on top face of T cell.
!-----------------------------------------------------------------------
 do k=1,nz-1
  do j=js_pe,je_pe
   do i=is_pe,ie_pe
!    eastward slopes at the top of T cells
     K_thkbi(i,j,k) = A_thkbi*cost(j)**(2*biha_friction_cosPower)
     do ip=0,1
      do kr=0,1 ! ddx2t(i-1+ip,j,k+kr,1)
        sxb = -drodx2b(i,j,k,ip,kr)/(min(0d0,drodzb(i,j,k,kr))-epsln)  ! i-1,k+1
        clip = taper(sxb,zw(k),mld(i,j)+mld0_thk)
        Ai_bx(i,j,k,ip,kr) =  clip*sxb*maskW(i,j,k)
        K_thkbi(i,j,k) = K_thkbi(i,j,k)*clip
      enddo
     enddo
!    northward slopes at the top of T cells
     do jp=0,1
      do kr=0,1
        syb = -drody2b(i,j,k,jp,kr)/(min(0d0,drodzb(i,j,k,kr))-epsln)
        clip = taper(syb,zw(k),mld(i,j)+mld0_thk)
        Ai_by(i,j,k,jp,kr) = clip*syb  *maskW(i,j,k)
        K_thkbi(i,j,k) = K_thkbi(i,j,k)*clip
      enddo
     enddo
   enddo
  enddo
 enddo
 Ai_bx(:,:,nz,:,:)=0d0
 Ai_by(:,:,nz,:,:)=0d0
 
 
 if (enable_conserve_energy) diss_biha_skew = 0d0
 
 do j=js_pe-onx,je_pe+onx
  diff_loc(:,j,:) = A_thkbi*cost(j)**(2*biha_friction_cosPower)
 enddo 
 
 call biharmonic_thickness_skew_divergence(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,&
                                           temp,Ai_ez,Ai_nz,Ai_bx,Ai_by,diff_loc,diss_biha_skew, .true.)
 call biharmonic_thickness_skew_divergence(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz, &
                                           salt,Ai_ez,Ai_nz,Ai_bx,Ai_by,diff_loc,diss_biha_skew, .false.)
  
  
 if (mod(itt,int(snapint/dt_tracer)) ==0 .or. itt==0) &
             call biharmonic_thickness_diag_store_tensor(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz, &
                                                         -Ai_by(:,:,:,0,0)*diff_loc,Ai_bx(:,:,:,0,0)*diff_loc)

  
 if (enable_biharmonic_thickness_backscatter_skew) &
   call biharmonic_thickness_backscatter_skew(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,& 
                                              drdTS,ddzt,ddxt,ddyt)
  
end subroutine biharmonic_thickness_skew_tensor




subroutine biharmonic_thickness_skew_divergence(is_,ie_,js_,je_,nz_,& 
                                                var,Ai_ez,Ai_nz,Ai_bx,Ai_by,diff_loc,diss_loc,istemp)
 use main_module 
 use biharmonic_thickness_module
 implicit none
 integer :: is_,ie_,js_,je_,nz_
 real*8 :: Ai_ez(is_:ie_,js_:je_,nz_,0:1,0:1) 
 real*8 :: Ai_nz(is_:ie_,js_:je_,nz_,0:1,0:1) 
 real*8 :: Ai_bx(is_:ie_,js_:je_,nz_,0:1,0:1)
 real*8 :: Ai_by(is_:ie_,js_:je_,nz_,0:1,0:1)
 real*8 :: var(is_:ie_,js_:je_,nz_,3)
 real*8 :: diss_loc(is_:ie_,js_:je_,nz_),diff_loc(is_:ie_,js_:je_,nz_)
 logical, intent(in) :: istemp
 integer :: i,j,k,kr,ip,jp,km1kr,kpkr,ks
 real*8 :: sumz,sumx,sumy,fxa
 real*8 :: aloc(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz)
 real*8 :: bloc(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz)

!-----------------------------------------------------------------------
!     construct total isoneutral tracer flux at east face of "T" cells 
!-----------------------------------------------------------------------
 do k=1,nz
  do j=js_pe,je_pe
   do i=is_pe-1,ie_pe
     !diffloc = A_thkbi*cost(j)**(2*biha_friction_cosPower)
     fxa = 0.5*(diff_loc(i,j,k) + diff_loc(i+1,j,k) )
     sumz = 0.
     do kr=0,1
       km1kr = max(k-1+kr,1)
       kpkr  = min(k+kr,nz)
       do ip=0,1
         sumz = sumz + fxa*Ai_ez(i,j,k,ip,kr) *(var(i+ip,j,kpkr,tau)-var(i+ip,j,km1kr,tau))
       enddo
     enddo
     flux_east(i,j,k) = sumz/(4*dzt(k))
   enddo
  enddo
 enddo
!-----------------------------------------------------------------------
!     construct total isoneutral tracer flux at north face of "T" cells 
!-----------------------------------------------------------------------
 do k=1,nz
  do j=js_pe-1,je_pe
   do i=is_pe,ie_pe
     !diffloc = A_thkbi*0.5*( cost(j)**(2*biha_friction_cosPower)+cost(j+1)**(2*biha_friction_cosPower))
     fxa = 0.5*(diff_loc(i,j,k) + diff_loc(i,j+1,k) )
     sumz    = 0.
     do kr=0,1
       km1kr = max(k-1+kr,1)
       kpkr  = min(k+kr,nz)
       do jp=0,1
         sumz = sumz + fxa*Ai_nz(i,j,k,jp,kr) *(var(i,j+jp,kpkr,tau)-var(i,j+jp,km1kr,tau))
       enddo
     enddo
     flux_north(i,j,k) = cosu(j)*sumz/(4*dzt(k))
   enddo
  enddo
 enddo
!-----------------------------------------------------------------------
!     compute the vertical tracer flux "flux_top" containing the K31
!     and K32 components which are to be solved explicitly. 
!----------------------------------------------------------------------- 
 do k=1,nz-1
  do j=js_pe,je_pe
   do i=is_pe,ie_pe
    !diffloc = A_thkbi*cost(j)**(2*biha_friction_cosPower)
    sumx = 0.
    do ip=0,1
     do kr=0,1
      sumx = sumx - diff_loc(i,j,k)*Ai_bx(i,j,k,ip,kr)/cost(j)*(var(i+ip,j,k+kr,tau) - var(i-1+ip,j,k+kr,tau))
     enddo
    enddo
    sumy    = 0.
    do jp=0,1
     do kr=0,1
      sumy = sumy - diff_loc(i,j,k)*Ai_by(i,j,k,jp,kr)*cosu(j-1+jp)* (var(i,j+jp,k+kr,tau)-var(i,j-1+jp,k+kr,tau))
     enddo
    enddo
    flux_top(i,j,k) = sumx/(4*dxt(i)) +sumy/(4*dyt(j)*cost(j) ) 
   enddo
  enddo
 enddo
 flux_top(:,:,nz)=0d0
  
 do j=js_pe,je_pe
   do i=is_pe,ie_pe
      aloc(i,j,:)=maskT(i,j,:)*( (flux_east(i,j,:)-  flux_east(i-1,j,:))/(cost(j)*dxt(i)) &
                                +(flux_north(i,j,:)- flux_north(i,j-1,:))/(cost(j)*dyt(j)) )
   enddo
 enddo
 
 k=1; aloc(:,:,k) = aloc(:,:,k) + maskT(:,:,k)*flux_top(:,:,k)/dzt(k)   
 do k=2,nz
   aloc(:,:,k) = aloc(:,:,k) + maskT(:,:,k)*(flux_top(:,:,k)- flux_top(:,:,k-1))/dzt(k)   
 enddo

 do j=js_pe,je_pe
   do i=is_pe,ie_pe
      var(i,j,:,taup1) = var(i,j,:,taup1) + dt_tracer*aloc(i,j,:)
   enddo
 enddo

!---------------------------------------------------------------------------------
! dissipation by isopycnal mixing
!---------------------------------------------------------------------------------
if (enable_conserve_energy) then

 if (istemp) then
  bloc(:,:,:) = int_drhodT(:,:,:,tau)
 else
  bloc(:,:,:) = int_drhodS(:,:,:,tau)
 endif

 do k=1,nz
   do j=js_pe-onx+1,je_pe+onx-1
    do i=is_pe-onx+1,ie_pe+onx-1
     fxa = bloc(i,j,k)
     aloc(i,j,k) =+0.5*grav/rho_0*( (bloc(i+1,j,k)-fxa)*flux_east(i  ,j,k) &
                                   +(fxa-bloc(i-1,j,k))*flux_east(i-1,j,k) ) /(dxt(i)*cost(j))  &
                  +0.5*grav/rho_0*( (bloc(i,j+1,k)-fxa)*flux_north(i,j  ,k) &
                                   +(fxa-bloc(i,j-1,k))*flux_north(i,j-1,k) ) /(dyt(j)*cost(j)) 
    enddo
   enddo
 end do
 !---------------------------------------------------------------------------------
 ! dissipation interpolated on W-grid
 !---------------------------------------------------------------------------------
 do j=js_pe-onx,je_pe+onx
   do i=is_pe-onx,ie_pe+onx
     ks=kbot(i,j)
     if (ks>0) then
      k=ks; diss_loc(i,j,k) = diss_loc(i,j,k)+ &
                         0.5*(aloc(i,j,k)+aloc(i,j,k+1)) + 0.5*aloc(i,j,k)*dzw(max(1,k-1))/dzw(k)
      do k=ks+1,nz-1
       diss_loc(i,j,k) = diss_loc(i,j,k)+ 0.5*(aloc(i,j,k) +aloc(i,j,k+1))
      enddo
      k=nz; diss_loc(i,j,k) = diss_loc(i,j,k)+ aloc(i,j,k)
     endif
   enddo
 enddo
 !---------------------------------------------------------------------------------
 ! dissipation by vertical component of skew mixing
 !---------------------------------------------------------------------------------
  do k=1,nz-1
   do j=js_pe,je_pe
    do i=is_pe,ie_pe
     fxa = (-bloc(i,j,k+1) +bloc(i,j,k))/dzw(k)
     diss_loc(i,j,k) = diss_loc(i,j,k)  -grav/rho_0*fxa*flux_top(i,j,k)*maskW(i,j,k)  
    enddo
   enddo
  end do
endif

end subroutine biharmonic_thickness_skew_divergence 

