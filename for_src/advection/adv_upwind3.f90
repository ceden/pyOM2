





subroutine adv_flux_upwind3(is_,ie_,js_,je_,nz_,adv_fe,adv_fn,adv_ft,var)
!---------------------------------------------------------------------------------
! upwind third order from MITgcm
! Calculates the area integrated zonal flux due to advection of a tracer
! using upwind biased third-order interpolation (or the $\kappa=1/3$ scheme):
! \begin{equation*}
! F^x_{adv} = U \overline{ \theta  - \frac{1}{6} \delta_{ii} \theta }^i
!                 + \frac{1}{12} |U| \delta_{iii} \theta
! \end{equation*}
! Near boundaries, mask all the gradients ==> still 3rd O.
!---------------------------------------------------------------------------------
      use main_module   
      implicit none
      integer, intent(in) :: is_,ie_,js_,je_,nz_
      real*8, intent(inout) :: adv_fe(is_:ie_,js_:je_,nz_), adv_fn(is_:ie_,js_:je_,nz_)
      real*8, intent(inout) :: adv_ft(is_:ie_,js_:je_,nz_),    var(is_:ie_,js_:je_,nz_)
      integer :: i,j,k,km1,kp2
      real*8 :: Rjp,Rj,Rjm,Rjjp,Rjjm
      real*8, parameter :: oneSixth=1.D0/6.D0

      do k=1,nz
       do j=js_pe,je_pe
        do i=is_pe-1,ie_pe
         Rjp=(var(i+2,j,k)-var(i+1,j,k))*maskU(i+1,j,k)
         Rj =(var(i+1,j,k)-var(i  ,j,k))*maskU(i  ,j,k)
         Rjm=(var(i  ,j,k)-var(i-1,j,k))*maskU(i-1,j,k)
         Rjjp=Rjp-Rj; Rjjm=Rj-Rjm
         adv_fe(i,j,k) = u(i,j,k,tau)*(var(i+1,j,k)+var(i,j,k) -oneSixth*( Rjjp+Rjjm )  )*0.5d0   &
                    +ABS(u(i,j,k,tau))*0.5d0*oneSixth*( Rjjp-Rjjm )
        enddo
       enddo
     enddo

      do k=1,nz
       do j=js_pe-1,je_pe
        do i=is_pe,ie_pe
         Rjp=(var(i,j+2,k)-var(i,j+1,k))*maskV(i,j+1,k)
         Rj =(var(i,j+1,k)-var(i,j  ,k))*maskV(i,j  ,k)
         Rjm=(var(i,j  ,k)-var(i,j-1,k))*maskV(i,j-1,k)
         Rjjp=Rjp-Rj; Rjjm=Rj-Rjm
         adv_fn(i,j,k) = cosu(j)*v(i,j,k,tau)*(var(i,j+1,k)+var(i,j,k)  -oneSixth*( Rjjp+Rjjm )  )*0.5d0   &
                    +ABS(cosu(j)*v(i,j,k,tau))*0.5d0*oneSixth*( Rjjp-Rjjm )
        enddo
       enddo
      enddo
 
      do k=1,nz-1
       kp2=min(nz,k+2); !if (kp2>np) kp2=3
       km1=max(1,k-1) !if (km1<1) km1=np-2
       do j=js_pe,je_pe
        do i=is_pe,ie_pe
         Rjp=(var(i,j,kp2)-var(i,j,k+1))*maskW(i,j,k+1)
         Rj =(var(i,j,k+1)-var(i,j,k  ))*maskW(i,j,k  )
         Rjm=(var(i,j,k  )-var(i,j,km1))*maskW(i,j,km1)
         Rjjp=Rjp-Rj; Rjjm=Rj-Rjm
         adv_ft(i,j,k) = w(i,j,k,tau)*(var(i,j,k+1)+var(i,j,k)   -oneSixth*( Rjjp+Rjjm )   )*0.5d0   &
                                +ABS(w(i,j,k,tau))*0.5d0*oneSixth*( Rjjp-Rjjm )
        enddo
       enddo
      enddo
      adv_ft(:,:,nz)=0.0
end subroutine adv_flux_upwind3

