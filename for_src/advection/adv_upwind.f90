

  

subroutine adv_flux_upwind_wgrid(is_,ie_,js_,je_,nz_,adv_fe,adv_fn,adv_ft,var)
!---------------------------------------------------------------------------------
! Calculates advection of a tracer defined on Wgrid
!---------------------------------------------------------------------------------
      use main_module   
      implicit none
      integer, intent(in) :: is_,ie_,js_,je_,nz_
      real*8, intent(inout) :: adv_fe(is_:ie_,js_:je_,nz_), adv_fn(is_:ie_,js_:je_,nz_)
      real*8, intent(inout) :: adv_ft(is_:ie_,js_:je_,nz_),    var(is_:ie_,js_:je_,nz_)
      integer :: i,j,k
      real*8 :: Rj
      real*8 :: maskUtr,maskVtr,maskWtr
      maskUtr(i,j,k) = maskW(i+1,j,k)*maskW(i,j,k)
      maskVtr(i,j,k) = maskW(i,j+1,k)*maskW(i,j,k)
      maskWtr(i,j,k) = maskW(i,j,k+1)*maskW(i,j,k)

      do k=1,nz
       do j=js_pe,je_pe
        do i=is_pe-1,ie_pe
         Rj =(var(i+1,j,k)-var(i  ,j,k))*maskUtr(i  ,j,k)
         adv_fe(i,j,k) = u_wgrid(i,j,k)*(var(i+1,j,k)+var(i,j,k))*0.5d0  -ABS(u_wgrid(i,j,k))*Rj*0.5d0
        enddo
       enddo
      enddo

      do k=1,nz
       do j=js_pe-1,je_pe
        do i=is_pe,ie_pe
         Rj =(var(i,j+1,k)-var(i,j  ,k))*maskVtr(i,j  ,k)
         adv_fn(i,j,k) = cosu(j)*v_wgrid(i,j,k)*(var(i,j+1,k)+var(i,j,k))*0.5d0 -ABS(cosu(j)*v_wgrid(i,j,k))*Rj*0.5d0
        enddo
       enddo
      enddo
 
      do k=1,nz-1
       do j=js_pe,je_pe
        do i=is_pe,ie_pe
         Rj =(var(i,j,k+1)-var(i,j,k  ))*maskWtr(i,j,k  )
         adv_ft(i,j,k) = w_wgrid(i,j,k)*(var(i,j,k+1)+var(i,j,k))*0.5d0 -ABS(w_wgrid(i,j,k))*Rj*0.5d0
        enddo
       enddo
      enddo
      adv_ft(:,:,nz)=0.0
end subroutine adv_flux_upwind_wgrid



