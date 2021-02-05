



subroutine momentum_advection
!=======================================================================
! Advection of momentum with second order which is energy conserving
!=======================================================================
  use main_module   
  implicit none
  integer :: i,j,k
  real*8 :: utr(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz)
  real*8 :: vtr(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz)
  real*8 :: wtr(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz)


 !---------------------------------------------------------------------------------
 !  Code from MITgcm
 !---------------------------------------------------------------------------------

!        uTrans(i,j) = u(i,j)*dyG(i,j)*drF(k)
!        vTrans(i,j) = v(i,j)*dxG(i,j)*drF(k)

!        fZon(i,j) = 0.25*( uTrans(i,j) + uTrans(i+1,j) ) *( u(i,j) + u(i+1,j) )
!        fMer(i,j) = 0.25*( vTrans(i,j) + vTrans(i-1,j) ) *( u(i,j) + u(i,j-1) )

!          gU(i,j,k,bi,bj) =  -
!     &     *( ( fZon(i,j  )  - fZon(i-1,j)  )
!     &       +( fMer(i,j+1)  - fMer(i,  j)  )
!     &       +( fVerUkp(i,j) - fVerUkm(i,j) )
!     &     ) /drF(k)   / rAw(i,j) 


!        fZon(i,j) = 0.25*( uTrans(i,j) + uTrans(i,j-1) )  *(v(i,j) + v(i-1,j) )
!        fMer(i,j) = 0.25*( vTrans(i,j) + vTrans(i,j+1) )  *(v(i,j) +  v(i,j+1) )

!          gV(i,j,k,bi,bj) =  -recip_drF(k)*recip_rAs(i,j,bi,bj)
!     &     *( ( fZon(i+1,j)  - fZon(i,j  )  )
!     &       +( fMer(i,  j)  - fMer(i,j-1)  )
!     &       +( fVerVkp(i,j) - fVerVkm(i,j) )
!     &     )


 do j=js_pe-onx,je_pe+onx
  do i=is_pe-onx,ie_pe+onx
    utr(i,j,:) = dzt(:)*dyt(j)*u(i,j,:,tau)*maskU(i,j,:)
    vtr(i,j,:) = dzt(:)*cosu(j)*dxt(i)*v(i,j,:,tau)*maskV(i,j,:)
    wtr(i,j,:) = area_t(i,j)*w(i,j,:,tau)*maskW(i,j,:)
  enddo
 enddo

 !---------------------------------------------------------------------------------
 ! for zonal momentum
 !---------------------------------------------------------------------------------
 do j=js_pe,je_pe
  do i=is_pe-1,ie_pe
    flux_east(i,j,:) = 0.25*(u(i,j,:,tau)+u(i+1,j,:,tau))*(utr(i+1,j,:)+utr(i,j,:))
  enddo
 enddo
 do j=js_pe-1,je_pe
  do i=is_pe,ie_pe
     flux_north(i,j,:) = 0.25*(u(i,j,:,tau)+u(i,j+1,:,tau))*(vtr(i+1,j,:)+vtr(i,j,:))
  enddo
 enddo
 do k=1,nz-1
  do j=js_pe,je_pe
   do i=is_pe,ie_pe
      flux_top(i,j,k) = 0.25*(u(i,j,k+1,tau)+u(i,j,k,tau))*(wtr(i,j,k)+wtr(i+1,j,k))
   enddo
  enddo
 enddo
 flux_top(:,:,nz)=0.0
 do j=js_pe,je_pe
  do i=is_pe,ie_pe
    du_adv(i,j,:) =  - maskU(i,j,:)*( flux_east(i,j,:) -flux_east(i-1,j,:) &
                                     +flux_north(i,j,:)-flux_north(i,j-1,:))/(area_u(i,j)*dzt(:))
  enddo
 enddo
 k=1; du_adv(:,:,k) = du_adv(:,:,k) - maskU(:,:,k)*flux_top(:,:,k)/(area_u(:,:)*dzt(k))
 do k=2,nz
     du_adv(:,:,k) = du_adv(:,:,k) - maskU(:,:,k)*(flux_top(:,:,k)-flux_top(:,:,k-1))/(area_u(:,:)*dzt(k))
 enddo
 !---------------------------------------------------------------------------------
 ! for meridional momentum
 !---------------------------------------------------------------------------------
 do j=js_pe,je_pe
  do i=is_pe-1,ie_pe
    flux_east(i,j,:) = 0.25*(v(i,j,:,tau)+v(i+1,j,:,tau))*(utr(i,j+1,:)+utr(i,j,:))
  enddo
 enddo
 do j=js_pe-1,je_pe
  do i=is_pe,ie_pe
     flux_north(i,j,:) = 0.25*(v(i,j,:,tau)+v(i,j+1,:,tau))*(vtr(i,j+1,:)+vtr(i,j,:))
  enddo
 enddo
 do k=1,nz-1
  do j=js_pe,je_pe
   do i=is_pe,ie_pe
      flux_top(i,j,k) = 0.25*(v(i,j,k+1,tau)+v(i,j,k,tau))*(wtr(i,j,k)+wtr(i,j+1,k))
   enddo
  enddo
 enddo
 flux_top(:,:,nz)=0.0
 do j=js_pe,je_pe
  do i=is_pe,ie_pe
    dv_adv(i,j,:) = - maskV(i,j,:)*( flux_east(i,j,:) -flux_east(i-1,j,:) &
                                    +flux_north(i,j,:)-flux_north(i,j-1,:))/(area_v(i,j)*dzt(:))
  enddo
 enddo
 k=1; dv_adv(:,:,k) = dv_adv(:,:,k) - maskV(:,:,k)*flux_top(:,:,k)/(area_v(:,:)*dzt(k))
 do k=2,nz
     dv_adv(:,:,k) = dv_adv(:,:,k) - maskV(:,:,k)*(flux_top(:,:,k)-flux_top(:,:,k-1))/(area_v(:,:)*dzt(k))
 enddo


 if (.not. enable_hydrostatic) then
  !---------------------------------------------------------------------------------
  ! for vertical momentum
  !---------------------------------------------------------------------------------
  do k=1,nz
   do j=js_pe,je_pe
    do i=is_pe-1,ie_pe
     flux_east(i,j,k) = 0.5*(w(i,j,k,tau)+w(i+1,j,k,tau))*(u(i,j,k,tau)+u(i,j,min(nz,k+1),tau))*0.5*maskW(i+1,j,k)*maskW(i,j,k)
    enddo
   enddo
  enddo
  do k=1,nz
   do j=js_pe-1,je_pe
    do i=is_pe,ie_pe
     flux_north(i,j,k) = 0.5*(w(i,j,k,tau)+w(i,j+1,k,tau))* & 
                         (v(i,j,k,tau)+v(i,j,min(nz,k+1),tau))*0.5*maskW(i,j+1,k)*maskW(i,j,k)*cosu(j)
    enddo
   enddo
  enddo
  do k=1,nz-1
   do j=js_pe,je_pe
    do i=is_pe,ie_pe
      flux_top(i,j,k) = 0.5*(w(i,j,k+1,tau)+w(i,j,k,tau))*(w(i,j,k,tau)+w(i,j,k+1,tau))*0.5*maskW(i,j,k+1)*maskW(i,j,k)
    enddo
   enddo
  enddo
  flux_top(:,:,nz)=0.0
  do j=js_pe,je_pe
   do i=is_pe,ie_pe
     dw_adv(i,j,:)=   maskW(i,j,:)* (-( flux_east(i,j,:)-  flux_east(i-1,j,:))/(cost(j)*dxt(i)) &
                                      -(flux_north(i,j,:)- flux_north(i,j-1,:))/(cost(j)*dyt(j)) )
   enddo
  enddo
  k=1; dw_adv(:,:,k) = dw_adv(:,:,k) - maskW(:,:,k)*flux_top(:,:,k)/dzw(k)
  do k=2,nz
     dw_adv(:,:,k) = dw_adv(:,:,k) - maskW(:,:,k)*(flux_top(:,:,k)-flux_top(:,:,k-1))/dzw(k)
  enddo
 endif

end subroutine momentum_advection



