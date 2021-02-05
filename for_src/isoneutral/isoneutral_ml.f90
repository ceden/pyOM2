


subroutine mixed_layer_para
!=======================================================================
! Parameterization of mixed layer eddies from Brüggemann and Eden (2014):
!  
! u'b' = - A H^2 f^3 alpha^2 L^2/f^2
! v'b' = - A H^2 f^3 alpha^2 M^2/f^2
! w'b' =   B H^2 f^3 alpha^2   
!  
! with alpha = Ro/delta = U/(fH), U=M^2 H/f -> alpha = sqrt(L^4+M^4)/f^2 
!  
! for FK: A = 2   C_F mu Ri         and  B = C_F mu  
! for ST: A = 8/5 C_S mu (1+Ri)^0.5 and  B = C_S mu 1/(1+Ri)^0.5 
!  
! streamfunction is given by
!  psi_x = - ( v'b' N^2 - w'b' M^2)/(N^4 + L^4 + M^4) 
!  psi_y = - (-u'b' N^2 + w'b' L^2)/(N^4 + L^4 + M^4)
!  psi_z = - ( u'b' M^2 - v'b' L^2)/(N^4 + L^4 + M^4)
!   
! Held + Schneider streamfct for N^2=0
!  psi_x =  B H^2 M^2/f
!  psi_y = -B H^2 L^2/f
!  psi_z =  0 
!  
! general streamfct.:
!  psi_x =  (H^2 M^2/f) (B+ A N^2/f^2) (L^4+M^4)/(N^4+L^4+M^4) 
!  psi_y = -(H^2 L^2/f) (B+ A N^2/f^2) (L^4+M^4)/(N^4+L^4+M^4) 
!  psi_z = 0
!  
! bolus velocity by
!  vn x psi =  (dx,dy,dz) x (psi_x,psi_y,psi_z) = (dy psi_z - dz psi_y , dz psi_x - dx psi_z, dx psi_y - dy psi_x)
!  u^* =  (psi_z)_y - (psi_y)_z 
!  v^* = -(psi_z)_x + (psi_x)_z
!  
! Diffusivity is given by
!   K = - (u'b' L^2 + v'b'M^2 + w'b'N^2)/(L^4+M^4+N^4)
!     =  H^2/f (  A L^4/f^2 + A M^4/f^2 - B N^2 )    (L^4+M^4)  /(L^4+M^4+N^4)
!=======================================================================
 use main_module   
 use isoneutral_module
 implicit none
 integer :: i,j,k
 real*8 :: mldu(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx)
 real*8 :: mldv(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx)
 
 real*8 :: Lsqr(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz)
 real*8 :: Msqr(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz)
 
 real*8 :: Nsqrvint(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx)
 real*8 :: Lsqrvint(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx)
 real*8 :: Msqrvint(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx)
 real*8 :: L4M4(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx)
 
 real*8 :: fCorTau(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx)
 real*8 :: Aloc(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz)
 real*8 :: Bloc(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz)
 real*8 :: mu(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx)
 real*8 :: fxa,get_rho,N2_min,N2_int
  !-----------------------------------------------------------------------
  ! find mixed layer depth 
  !-----------------------------------------------------------------------
  if (ml_algorithm == 1) then
    ! as given in Fox-Kemper et al. 
    mld=dzw(nz)/2.
    do j=js_pe-onx,je_pe+onx
     do i=is_pe-onx,ie_pe+onx
       if (kbot(i,j)>0) then
         N2_min = 1e18; N2_int=0d0; 
         do k=nz-1,kbot(i,j),-1
            fxa = max(0d0,Nsqr(i,j,k,tau))
            if ( fxa-N2_min  .gt. ml_Nsqr_fac*N2_int/mld(i,j)) exit
            N2_min = min(fxa,N2_min)
            N2_int = N2_int + fxa*dzw(k)
            mld(i,j)  = mld(i,j) + dzw(k)
         enddo
       endif
     enddo
    enddo
 elseif (ml_algorithm == 3) then
    do j=js_pe-onx,je_pe+onx
     do i=is_pe-onx,ie_pe+onx
         if (kbot(i,j)>0) then
          do k=nz-1,kbot(i,j),-1
            mld(i,j) = zw(nz)-zw(k)
            if ( get_rho(salt(i,j,k,tau),temp(i,j,k,tau),abs(zt(nz)))  .gt.  rho(i,j,nz,tau) +ml_Delta_rho ) exit
          enddo
         endif
     enddo
    enddo
  else
    call halt_stop('ERROR in mixed_layer_para: wrong parameter for ml_algorithm ')
  endif

  mldu=mld; mldv=mld
  do i=is_pe-onx,ie_pe+onx-1
    mldu(i,:) = 0.5*(mld(i+1,:)+mld(i,:)) 
  enddo
  do j=js_pe-onx,je_pe+onx-1
    mldv(:,j) = 0.5*(mld(:,j+1)+mld(:,j)) 
  enddo
  !-----------------------------------------------------------------------
  ! calculate horizontal and vertical buoyancy gradients
  !-----------------------------------------------------------------------
  Lsqr=0D0; Msqr=0D0
  fxa =  -grav/rho_0
  do j=js_pe-onx,je_pe+onx
    do i=is_pe-onx,ie_pe+onx-1
      Lsqr(i,j,:) = fxa*(rho(i+1,j,:,tau)-rho(i,j,:,tau)) / (dxu(i)*cost(j)) * maskU(i,j,:)
    enddo
  enddo
  do j=js_pe-onx,je_pe+onx-1
    do i=is_pe-onx,ie_pe+onx
      Msqr(i,j,:) = fxa*(rho(i,j+1,:,tau)-rho(i,j,:,tau)) / dyu(j) * maskV(i,j,:)
    enddo
  enddo
  fCorTau = sqrt( coriolis_t**2 + ml_invTau**2 )
  !-----------------------------------------------------------------------
  ! integrate buoyancy gradients over mixed layer depth 
  !-----------------------------------------------------------------------
  Lsqrvint = 0.0
  do k=1,nz
     Lsqrvint = Lsqrvint + Lsqr(:,:,k)*  max(0d0, min(dzt(k),mldu+zw(k)))
  enddo
  Msqrvint = 0.0
  do k=1,nz
    Msqrvint = Msqrvint + Msqr(:,:,k)* max(0D0, min(dzt(k), mldv+zw(k)))
  enddo
  call border_exchg_xy(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,Lsqrvint)
  call setcyclic_xy   (is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,Lsqrvint)
  call border_exchg_xy(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,Msqrvint)
  call setcyclic_xy   (is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,Msqrvint)
  Nsqrvint = 0.0
  do k=1,nz-1
    do j=js_pe-onx,je_pe+onx
      do i=is_pe-onx,ie_pe+onx
        fxa = 0.0
        if ( mld(i,j)+zt(k+1) .ge. 0.0 ) fxa = dzw(k)
        Nsqrvint(i,j) = Nsqrvint(i,j) + max(Nsqr(i,j,k,tau),0d0)*fxa* maskW(i,j,k) 
      enddo
    enddo
  enddo
  !-----------------------------------------------------------------------
  ! calculate factor (L^4 + M^4)/(N^4+M^4+L^4) on T grid
  !-----------------------------------------------------------------------
  L4M4=0.0
  do j=js_pe,je_pe
    do i=is_pe,ie_pe
       L4M4(i,j) = (0.5*(Msqrvint(i,j-1)+Msqrvint(i,j)))**2 +  (0.5*(Lsqrvint(i-1,j)+Lsqrvint(i,j)))**2
       L4M4(i,j) = L4M4(i,j)/(Nsqrvint(i,j)**2 + L4M4(i,j) +1d-18 )
    enddo
  enddo
  call border_exchg_xy(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,L4M4)
  call setcyclic_xy   (is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,L4M4)

  !-----------------------------------------------------------------------
  ! compute Ri   
  !-----------------------------------------------------------------------
  do j=js_pe-onx+1,je_pe+onx
    do i=is_pe-onx+1,ie_pe+onx
      Ri_ml(i,j) = Nsqrvint(i,j)*coriolis_t(i,j)**2 * mld(i,j) & 
                * 1.0 / ( 0.25*(Lsqrvint(i-1,j)+Lsqrvint(i,j))**2 &
                        + 0.25*(Msqrvint(i,j-1)+Msqrvint(i,j))**2 + 1d-12)
    enddo
  enddo
  call border_exchg_xy(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,Ri_ml)
  call setcyclic_xy   (is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,Ri_ml)
 
  !-----------------------------------------------------------------------
  ! calculate factor B: for Fox-Kemper  B= C_F mu
  !                     for Stone       B= C_S mu /(1+Ri)^0.5 
  !-----------------------------------------------------------------------
  do k=1,nz
     mu = 0.0
     where ( abs(zt(k)) .le. mld ) 
        mu = zw(max(1,k-1))/(mld+dzw(k))
        mu = -4.0*mu*(mu+1)
     end where
     Bloc(:,:,k)=ml_Ce*mu
     if ( enable_ml_para_stone ) Bloc(:,:,k)=Bloc(:,:,k)/sqrt(1.0+Ri_ml)
  enddo

  !-----------------------------------------------------------------------
  ! calculate factor A: for Fox-Kemper  A = 2   C_F mu Ri 
  !                     for Stone       A = 8/5 C_S mu (1+Ri)^0.5  
  !-----------------------------------------------------------------------
  if (.not. enable_ml_HS_streamfct) then
    do k=1,nz
       mu = 0.0
       where ( abs(zt(k)) .le. mld ) 
          mu = zw(max(1,k-1))/(mld+dzw(k))
          mu = -4.0*mu*(mu+1)
       end where
       Aloc(:,:,k)=ml_Ce*mu

       if ( enable_ml_para_stone ) then
          Aloc(:,:,k)=8./5.*Aloc(:,:,k)*sqrt(1+Ri_ml)
       else
          Aloc(:,:,k)=2*Aloc(:,:,k)*Ri_ml
       endif
    enddo
  endif

  if (enable_ml_HS_streamfct) then
  !-----------------------------------------------------------------------
  !    compute simple Held+Schneider version of ML eddy streamfunction
  !    psi_y = -B H^2 L^2/f  and psi_x =  B H^2 M^2/f 
  !-----------------------------------------------------------------------
   do k=1,nz
    do j=js_pe-1,je_pe
      do i=is_pe-1,ie_pe
       B2_ml(i,j,k) = -(Bloc(i,j,k)+Bloc(i+1,j,k))*Lsqrvint(i,j)*mldu(i,j)/( fCorTau(i+1,j)+fCorTau(i,j) )
       B1_ml(i,j,k) =  (Bloc(i,j,k)+Bloc(i,j+1,k))*Msqrvint(i,j)*mldv(i,j)/( fCorTau(i,j+1)+fCorTau(i,j) )
      enddo
    enddo
  enddo
 else
  !-----------------------------------------------------------------------
  ! compute general version of ML eddy streamfunction
  !      psi_y = -(H^2 L^2/f) (B+ A N^2/f^2) (L^4+M^4)/(N^4+L^4+M^4) 
  !      psi_x =  (H^2 M^2/f) (B+ A N^2/f^2) (L^4+M^4)/(N^4+L^4+M^4) 
  !-----------------------------------------------------------------------
   do k=1,nz
    do j=js_pe-1,je_pe
      do i=is_pe-1,ie_pe
        fxa = (Nsqrvint(i+1,j)+Nsqrvint(i,j))/max(1d-6,mldu(i,j))/(fCorTau(i+1,j)+fCorTau(i,j))**2
        B2_ml(i,j,k) = -Lsqrvint(i,j)*mldu(i,j)/( fCorTau(i+1,j)+fCorTau(i,j) ) &
                       *(Bloc(i,j,k)+Bloc(i+1,j,k)+(Aloc(i,j,k)+Aloc(i+1,j,k))*fxa)*0.5*(L4M4(i+1,j)+L4M4(i,j))
        
        fxa = (Nsqrvint(i,j+1)+Nsqrvint(i,j))/max(1d-6,mldv(i,j))/(fCorTau(i,j+1)+fCorTau(i,j))**2
        B1_ml(i,j,k) = Msqrvint(i,j)*mldv(i,j)/( fCorTau(i,j+1)+fCorTau(i,j) ) &
                       *(Bloc(i,j,k)+Bloc(i,j+1,k)+(Aloc(i,j,k)+Aloc(i,j+1,k))*fxa)*0.5*(L4M4(i,j+1)+L4M4(i,j))

     enddo
    enddo
  enddo
  !-----------------------------------------------------------------------
  ! calculate also diffusivity
  !      K =  (  A H^2 L^4/f^3 + A H^2 M^4/f^3 - B H^2 N^2/f )  (L^4+M^4) /(L^4+M^4+N^4)
  ! very small, only for diagnostics
  !-----------------------------------------------------------------------
   do k=1,nz
    do j=js_pe,je_pe
      do i=is_pe,ie_pe
        K_ml(i,j,k)  = (0.5*(Lsqrvint(i,j)+Lsqrvint(i-1,j)))**2 +  (0.5*(Msqrvint(i,j)+Msqrvint(i,j-1)))**2
        K_ml(i,j,k)  = Aloc(i,j,k)*K_ml(i,j,k)/fCorTau(i,j)**3 - Bloc(i,j,k)*Nsqrvint(i,j)*mld(i,j)/fCorTau(i,j)
        K_ml(i,j,k)  = K_ml(i,j,k)*L4M4(i,j)
     enddo
    enddo
  enddo

 endif

  !-----------------------------------------------------------------------
  ! calculate eddy velocities from streamfct.
  !-----------------------------------------------------------------------

  do k=1,nz-1
    u_ml(:,:,k) = -(B2_ml(:,:,k+1)-B2_ml(:,:,k))/dzw(k)*maskU(:,:,k)
    v_ml(:,:,k) =  (B1_ml(:,:,k+1)-B1_ml(:,:,k))/dzw(k)*maskV(:,:,k)
  enddo
  k=nz
  u_ml(:,:,k) = -(-B2_ml(:,:,k))/dzw(k)*maskU(:,:,k)
  v_ml(:,:,k) =  (-B1_ml(:,:,k))/dzw(k)*maskV(:,:,k)
  !-----------------------------------------------------------------------
  ! calculate vertical eddy velocity by continuity equation
  !-----------------------------------------------------------------------
   w_ml=0d0
   k=1
   do j=js_pe,je_pe
     do i=is_pe,ie_pe
        w_ml(i,j,k) =-maskW(i,j,k)*dzt(k)* &
               ((        u_ml(i,j,k)-          u_ml(i-1,j,k))/(cost(j)*dxt(i)) &
               +(cosu(j)*v_ml(i,j,k)-cosu(j-1)*v_ml(i,j-1,k))/(cost(j)*dyt(j)) )
     enddo
   enddo
   do k=2,nz
     do j=js_pe,je_pe
      do i=is_pe,ie_pe
          w_ml(i,j,k) = w_ml(i,j,k-1)-maskW(i,j,k)*dzt(k)* &
               ((        u_ml(i,j,k)          -u_ml(i-1,j,k))/(cost(j)*dxt(i)) &
               +(cosu(j)*v_ml(i,j,k)-cosu(j-1)*v_ml(i,j-1,k))/(cost(j)*dyt(j)) )
      enddo
     enddo
   enddo

  !-----------------------------------------------------------------------
  ! calculate advection 
  !-----------------------------------------------------------------------
   call advect_with_bolus_flow(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,temp(:,:,:,tau),dtemp_ml,u_ml,v_ml,w_ml)
   call advect_with_bolus_flow(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,salt(:,:,:,tau),dsalt_ml,u_ml,v_ml,w_ml)
   do j=js_pe,je_pe
    do i=is_pe,ie_pe
      temp(i,j,:,taup1) = temp(i,j,:,taup1) + dt_tracer*dtemp_ml(i,j,:)
      salt(i,j,:,taup1) = salt(i,j,:,taup1) + dt_tracer*dsalt_ml(i,j,:)
     enddo
   enddo

end subroutine mixed_layer_para


subroutine mixed_layer_para_trac(is_,ie_,js_,je_,nz_,tr)
!-----------------------------------------------------------------------
! apply mixed layer parameterisation also to tracer
!-----------------------------------------------------------------------
 use main_module   
 use isoneutral_module   
 implicit none
 integer :: is_,ie_,js_,je_,nz_
 real*8 :: tr(is_:ie_,js_:je_,nz_,3)
 real*8 :: dtr(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz)
 integer :: i,j
 
 call advect_with_bolus_flow(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,tr(:,:,:,tau),dtr,u_ml,v_ml,w_ml)
 do j=js_pe,je_pe
  do i=is_pe,ie_pe
    tr(i,j,:,taup1) = tr(i,j,:,taup1) + dt_tracer*dtr(i,j,:)
   enddo
 enddo
end subroutine mixed_layer_para_trac





subroutine advect_with_bolus_flow(is_,ie_,js_,je_,nz_,tr,dtr,uadv,vadv,wadv)
!-----------------------------------------------------------------------
! calculate 2nd order advection tendency with bolus flow 
!-----------------------------------------------------------------------
  use main_module   
  implicit none
  integer :: is_,ie_,js_,je_,nz_
  real*8 :: uadv(is_:ie_,js_:je_,nz_), vadv(is_:ie_,js_:je_,nz_), wadv(is_:ie_,js_:je_,nz_)
  real*8 :: tr(is_:ie_,js_:je_,nz_),  dtr(is_:ie_,js_:je_,nz_)
  integer :: i,j,k

  do k=1,nz
   do j=js_pe,je_pe
    do i=is_pe-1,ie_pe
     flux_east(i,j,k)=0.5*(tr(i,j,k) + tr(i+1,j,k) )*uadv(i,j,k)*maskU(i,j,k)
    enddo
   enddo
  enddo
  do k=1,nz
   do j=js_pe-1,je_pe
    do i=is_pe,ie_pe
     flux_north(i,j,k)=cosu(j)*0.5*( tr(i,j,k) + tr(i,j+1,k) )*vadv(i,j,k)*maskV(i,j,k)
    enddo
   enddo
  enddo
  do k=1,nz-1
   do j=js_pe,je_pe
    do i=is_pe,ie_pe
     flux_top(i,j,k)=0.5*( tr(i,j,k) + tr(i,j,k+1) )*wadv(i,j,k)*maskW(i,j,k)
    enddo
   enddo
  enddo
  flux_top(:,:,nz)=0.0

  do j=js_pe,je_pe
    do i=is_pe,ie_pe
       dtr(i,j,:)=maskT(i,j,:)* (-( flux_east(i,j,:)-  flux_east(i-1,j,:))/(cost(j)*dxt(i)) &
                                 -( flux_north(i,j,:)- flux_north(i,j-1,:))/(cost(j)*dyt(j)) )
    enddo
  enddo
  k=1; dtr(:,:,k)=dtr(:,:,k)-maskT(:,:,k)*flux_top(:,:,k)/dzt(k)   
  do k=2,nz
    dtr(:,:,k)=dtr(:,:,k)-maskT(:,:,k)*(flux_top(:,:,k)- flux_top(:,:,k-1))/dzt(k)   
  enddo
end subroutine advect_with_bolus_flow
