


subroutine adv_flux_dst3(is_,ie_,js_,je_,nz_,adv_fe,adv_fn,adv_ft,var)
!---------------------------------------------------------------------------------
! from MITgcm
!   Compute advective Flux of Tracer using 3rd Order DST Sceheme with flux limiting               
!---------------------------------------------------------------------------------
      use main_module   
      implicit none
      integer, intent(in) :: is_,ie_,js_,je_,nz_
      real*8, intent(inout) :: adv_fe(is_:ie_,js_:je_,nz_), adv_fn(is_:ie_,js_:je_,nz_)
      real*8, intent(inout) :: adv_ft(is_:ie_,js_:je_,nz_),    var(is_:ie_,js_:je_,nz_)
      integer :: i,j,k,km1,kp2
      real*8 :: Rjp,Rj,Rjm,uCFL,d0,d1,thetaP,thetaM,psiP,psiM
      real*8, parameter :: oneSixth=1.D0/6.D0
      real*8, parameter :: thetaMax = 1.D+20 

      do k=1,nz
       do j=js_pe,je_pe
        do i=is_pe-1,ie_pe
         uCFL = ABS( u(i,j,k,tau)*dt_tracer/(cost(j)*dxt(i)) )
         Rjp=(var(i+2,j,k)-var(i+1,j,k))*maskU(i+1,j,k)
         Rj =(var(i+1,j,k)-var(i  ,j,k))*maskU(i  ,j,k)
         Rjm=(var(i  ,j,k)-var(i-1,j,k))*maskU(i-1,j,k)

         d0=(2d0 -uCFL)*(1d0 -uCFL)*oneSixth
         d1=(1d0 -uCFL*uCFL)*oneSixth
         IF ( ABS(Rj)*thetaMax .LE. ABS(Rjm) ) THEN
             thetaP=SIGN(thetaMax,Rjm*Rj)
         ELSE
          thetaP=Rjm/Rj
         ENDIF
         IF ( ABS(Rj)*thetaMax .LE. ABS(Rjp) ) THEN
          thetaM=SIGN(thetaMax,Rjp*Rj)
         ELSE
          thetaM=Rjp/Rj
         ENDIF
         psiP=d0+d1*thetaP
         psiP=MAX(0d0,MIN(MIN(1d0,psiP), thetaP*(1d0 -uCFL)/(uCFL+1d-20) ))
         psiM=d0+d1*thetaM
         psiM=MAX(0d0,MIN(MIN(1d0,psiM), thetaM*(1d0 -uCFL)/(uCFL+1d-20) ))

         adv_fe(i,j,k)= 0.5*(u(i,j,k,tau)+ABS(u(i,j,k,tau))) *( var(i  ,j,k) + psiP*Rj )  &
                       +0.5*(u(i,j,k,tau)-ABS(u(i,j,k,tau))) *( var(i+1,j,k) - psiM*Rj )
        enddo
       enddo
    enddo

      do k=1,nz
       do j=js_pe-1,je_pe
        do i=is_pe,ie_pe
         Rjp=(var(i,j+2,k)-var(i,j+1,k))*maskV(i,j+1,k)
         Rj =(var(i,j+1,k)-var(i,j  ,k))*maskV(i,j  ,k)
         Rjm=(var(i,j  ,k)-var(i,j-1,k))*maskV(i,j-1,k)
         uCFL = ABS( cosu(j)*v(i,j,k,tau)*dt_tracer/(cost(j)*dyt(j)) )

         d0=(2d0 -uCFL)*(1d0 -uCFL)*oneSixth
         d1=(1d0 -uCFL*uCFL)*oneSixth
         IF ( ABS(Rj)*thetaMax .LE. ABS(Rjm) ) THEN
             thetaP=SIGN(thetaMax,Rjm*Rj)
         ELSE
          thetaP=Rjm/Rj
         ENDIF
         IF ( ABS(Rj)*thetaMax .LE. ABS(Rjp) ) THEN
          thetaM=SIGN(thetaMax,Rjp*Rj)
         ELSE
          thetaM=Rjp/Rj
         ENDIF
         psiP=d0+d1*thetaP
         psiP=MAX(0d0,MIN(MIN(1d0,psiP), thetaP*(1d0 -uCFL)/(uCFL+1d-20) ))
         psiM=d0+d1*thetaM
         psiM=MAX(0d0,MIN(MIN(1d0,psiM), thetaM*(1d0 -uCFL)/(uCFL+1d-20) ))
         
         adv_fn(i,j,k)= 0.5*(v(i,j,k,tau)+ABS(v(i,j,k,tau))) *( var(i,j  ,k) + psiP*Rj )  &
                       +0.5*(v(i,j,k,tau)-ABS(v(i,j,k,tau))) *( var(i,j+1,k) - psiM*Rj )
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
         uCFL = ABS( w(i,j,k,tau)*dt_tracer/dzt(k) )

         d0=(2d0 -uCFL)*(1d0 -uCFL)*oneSixth
         d1=(1d0 -uCFL*uCFL)*oneSixth
         IF ( ABS(Rj)*thetaMax .LE. ABS(Rjm) ) THEN
             thetaP=SIGN(thetaMax,Rjm*Rj)
         ELSE
          thetaP=Rjm/Rj
         ENDIF
         IF ( ABS(Rj)*thetaMax .LE. ABS(Rjp) ) THEN
          thetaM=SIGN(thetaMax,Rjp*Rj)
         ELSE
          thetaM=Rjp/Rj
         ENDIF
         psiP=d0+d1*thetaP
         psiP=MAX(0d0,MIN(MIN(1d0,psiP), thetaP*(1d0 -uCFL)/(uCFL+1d-20) ))
         psiM=d0+d1*thetaM
         psiM=MAX(0d0,MIN(MIN(1d0,psiM), thetaM*(1d0 -uCFL)/(uCFL+1d-20) ))
         
         adv_ft(i,j,k)= 0.5*(w(i,j,k,tau)+ABS(w(i,j,k,tau))) *( var(i,j,k  ) + psiP*Rj )  &
                       +0.5*(w(i,j,k,tau)-ABS(w(i,j,k,tau))) *( var(i,j,k+1) - psiM*Rj )
        enddo
       enddo
      enddo
      adv_ft(:,:,nz)=0.0
end subroutine adv_flux_dst3








subroutine adv_flux_dst3_wgrid(is_,ie_,js_,je_,nz_,adv_fe,adv_fn,adv_ft,var)
!---------------------------------------------------------------------------------
! from MITgcm
!   Compute advective Flux of Tracer using 3rd Order DST Sceheme with flux limiting               
!---------------------------------------------------------------------------------
      use main_module   
      implicit none
      integer, intent(in) :: is_,ie_,js_,je_,nz_
      real*8, intent(inout) :: adv_fe(is_:ie_,js_:je_,nz_), adv_fn(is_:ie_,js_:je_,nz_)
      real*8, intent(inout) :: adv_ft(is_:ie_,js_:je_,nz_),    var(is_:ie_,js_:je_,nz_)
      integer :: i,j,k,km1,kp2
      real*8 :: Rjp,Rj,Rjm,uCFL,d0,d1,thetaP,thetaM,psiP,psiM
      real*8, parameter :: oneSixth=1.D0/6.D0
      real*8, parameter :: thetaMax = 1.D+20 
      real*8 :: maskUtr,maskVtr,maskWtr
      maskUtr(i,j,k) = maskW(i+1,j,k)*maskW(i,j,k)
      maskVtr(i,j,k) = maskW(i,j+1,k)*maskW(i,j,k)
      maskWtr(i,j,k) = maskW(i,j,k+1)*maskW(i,j,k)

      do k=1,nz
       do j=js_pe,je_pe
        do i=is_pe-1,ie_pe
         uCFL = ABS( u_wgrid(i,j,k)*dt_tracer/(cost(j)*dxt(i)) )
         Rjp=(var(i+2,j,k)-var(i+1,j,k))*maskUtr(i+1,j,k)
         Rj =(var(i+1,j,k)-var(i  ,j,k))*maskUtr(i  ,j,k)
         Rjm=(var(i  ,j,k)-var(i-1,j,k))*maskUtr(i-1,j,k)

         d0=(2d0 -uCFL)*(1d0 -uCFL)*oneSixth
         d1=(1d0 -uCFL*uCFL)*oneSixth
         IF ( ABS(Rj)*thetaMax .LE. ABS(Rjm) ) THEN
             thetaP=SIGN(thetaMax,Rjm*Rj)
         ELSE
          thetaP=Rjm/Rj
         ENDIF
         IF ( ABS(Rj)*thetaMax .LE. ABS(Rjp) ) THEN
          thetaM=SIGN(thetaMax,Rjp*Rj)
         ELSE
          thetaM=Rjp/Rj
         ENDIF
         psiP=d0+d1*thetaP
         psiP=MAX(0d0,MIN(MIN(1d0,psiP), thetaP*(1d0 -uCFL)/(uCFL+1d-20) ))
         psiM=d0+d1*thetaM
         psiM=MAX(0d0,MIN(MIN(1d0,psiM), thetaM*(1d0 -uCFL)/(uCFL+1d-20) ))

         adv_fe(i,j,k)= 0.5*(u_wgrid(i,j,k)+ABS(u_wgrid(i,j,k))) *( var(i  ,j,k) + psiP*Rj )  &
                       +0.5*(u_wgrid(i,j,k)-ABS(u_wgrid(i,j,k))) *( var(i+1,j,k) - psiM*Rj )
        enddo
       enddo
    enddo

      do k=1,nz
       do j=js_pe-1,je_pe
        do i=is_pe,ie_pe
         Rjp=(var(i,j+2,k)-var(i,j+1,k))*maskVtr(i,j+1,k)
         Rj =(var(i,j+1,k)-var(i,j  ,k))*maskVtr(i,j  ,k)
         Rjm=(var(i,j  ,k)-var(i,j-1,k))*maskVtr(i,j-1,k)
         uCFL = ABS( cosu(j)*v_wgrid(i,j,k)*dt_tracer/(cost(j)*dyt(j)) )

         d0=(2d0 -uCFL)*(1d0 -uCFL)*oneSixth
         d1=(1d0 -uCFL*uCFL)*oneSixth
         IF ( ABS(Rj)*thetaMax .LE. ABS(Rjm) ) THEN
             thetaP=SIGN(thetaMax,Rjm*Rj)
         ELSE
          thetaP=Rjm/Rj
         ENDIF
         IF ( ABS(Rj)*thetaMax .LE. ABS(Rjp) ) THEN
          thetaM=SIGN(thetaMax,Rjp*Rj)
         ELSE
          thetaM=Rjp/Rj
         ENDIF
         psiP=d0+d1*thetaP
         psiP=MAX(0d0,MIN(MIN(1d0,psiP), thetaP*(1d0 -uCFL)/(uCFL+1d-20) ))
         psiM=d0+d1*thetaM
         psiM=MAX(0d0,MIN(MIN(1d0,psiM), thetaM*(1d0 -uCFL)/(uCFL+1d-20) ))
         
         adv_fn(i,j,k)= 0.5*(v_wgrid(i,j,k)+ABS(v_wgrid(i,j,k))) *( var(i,j  ,k) + psiP*Rj )  &
                       +0.5*(v_wgrid(i,j,k)-ABS(v_wgrid(i,j,k))) *( var(i,j+1,k) - psiM*Rj )
        enddo
       enddo
      enddo
 
      do k=1,nz-1
       kp2=min(nz,k+2); !if (kp2>np) kp2=3
       km1=max(1,k-1) !if (km1<1) km1=np-2
       do j=js_pe,je_pe
        do i=is_pe,ie_pe
         Rjp=(var(i,j,kp2)-var(i,j,k+1))*maskWtr(i,j,k+1)
         Rj =(var(i,j,k+1)-var(i,j,k  ))*maskWtr(i,j,k  )
         Rjm=(var(i,j,k  )-var(i,j,km1))*maskWtr(i,j,km1)
         uCFL = ABS( w_wgrid(i,j,k)*dt_tracer/dzw(k) )

         d0=(2d0 -uCFL)*(1d0 -uCFL)*oneSixth
         d1=(1d0 -uCFL*uCFL)*oneSixth
         IF ( ABS(Rj)*thetaMax .LE. ABS(Rjm) ) THEN
             thetaP=SIGN(thetaMax,Rjm*Rj)
         ELSE
          thetaP=Rjm/Rj
         ENDIF
         IF ( ABS(Rj)*thetaMax .LE. ABS(Rjp) ) THEN
          thetaM=SIGN(thetaMax,Rjp*Rj)
         ELSE
          thetaM=Rjp/Rj
         ENDIF
         psiP=d0+d1*thetaP
         psiP=MAX(0d0,MIN(MIN(1d0,psiP), thetaP*(1d0 -uCFL)/(uCFL+1d-20) ))
         psiM=d0+d1*thetaM
         psiM=MAX(0d0,MIN(MIN(1d0,psiM), thetaM*(1d0 -uCFL)/(uCFL+1d-20) ))
         
         adv_ft(i,j,k)= 0.5*(w_wgrid(i,j,k)+ABS(w_wgrid(i,j,k))) *( var(i,j,k  ) + psiP*Rj )  &
                       +0.5*(w_wgrid(i,j,k)-ABS(w_wgrid(i,j,k))) *( var(i,j,k+1) - psiM*Rj )
        enddo
       enddo
      enddo
      adv_ft(:,:,nz)=0.0
end subroutine adv_flux_dst3_wgrid







