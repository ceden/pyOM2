

  

subroutine adv_flux_2nd(is_,ie_,js_,je_,nz_,adv_fe,adv_fn,adv_ft,var)
!---------------------------------------------------------------------------------
!      2nd order advective tracer flux
!---------------------------------------------------------------------------------
      use main_module   
      implicit none
      integer, intent(in) :: is_,ie_,js_,je_,nz_
      real*8, intent(inout) :: adv_fe(is_:ie_,js_:je_,nz_), adv_fn(is_:ie_,js_:je_,nz_)
      real*8, intent(inout) :: adv_ft(is_:ie_,js_:je_,nz_),    var(is_:ie_,js_:je_,nz_)
      integer :: i,j,k

      do k=1,nz
       do j=js_pe,je_pe
        do i=is_pe-1,ie_pe
         adv_fe(i,j,k)=0.5*(var(i,j,k) + var(i+1,j,k) )*u(i,j,k,tau)*maskU(i,j,k)
        enddo
       enddo
      enddo
      do k=1,nz
       do j=js_pe-1,je_pe
        do i=is_pe,ie_pe
         adv_fn(i,j,k)=cosu(j)*0.5*( var(i,j,k) + var(i,j+1,k) )*v(i,j,k,tau)*maskV(i,j,k)
        enddo
       enddo
      enddo
      do k=1,nz-1
       do j=js_pe,je_pe
        do i=is_pe,ie_pe
         adv_ft(i,j,k)=0.5*( var(i,j,k) + var(i,j,k+1) )*w(i,j,k,tau)*maskW(i,j,k)
        enddo
       enddo
      enddo
      adv_ft(:,:,nz)=0.0
end subroutine adv_flux_2nd



