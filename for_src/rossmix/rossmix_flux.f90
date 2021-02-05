



subroutine adv_flux_rossmix(is_,ie_,js_,je_,nphi_,adv_fe,adv_fn,adv_ft,var,nn,dtloc,uvel,vvel,wvel)
     implicit none
      integer, intent(in) :: is_,ie_,js_,je_,nphi_,nn
      real*8, intent(inout) :: adv_fe(is_:ie_,js_:je_,nphi_), adv_fn(is_:ie_,js_:je_,nphi_)
      real*8, intent(inout) :: adv_ft(is_:ie_,js_:je_,nphi_),    var(is_:ie_,js_:je_,nphi_,3)
      real*8, intent(in) :: uvel(is_:ie_,js_:je_,nphi_),vvel(is_:ie_,js_:je_,nphi_),wvel(is_:ie_,js_:je_,nphi_)
      real*8 :: dtloc
      call adv_flux_rossmix_superbee(is_,ie_,js_,je_,nphi_,adv_fe,adv_fn,adv_ft,var,nn,dtloc,uvel,vvel,wvel)
      !call adv_flux_rossmix_upwind(is_,ie_,js_,je_,nphi_,adv_fe,adv_fn,adv_ft,var,nn,dtloc,uvel,vvel,wvel)
end subroutine adv_flux_rossmix


subroutine adv_flux_rossmix_superbee(is_,ie_,js_,je_,nphi_,adv_fe,adv_fn,adv_ft,var,nn,ddt,uvel,vvel,wvel)
!=======================================================================
! Calculates advection of energy in physical and phase space
!=======================================================================
      use main_module   
      use rossmix_module   
      implicit none
      integer, intent(in) :: is_,ie_,js_,je_,nphi_,nn
      real*8, intent(inout) :: adv_fe(is_:ie_,js_:je_,nphi_), adv_fn(is_:ie_,js_:je_,nphi_)
      real*8, intent(inout) :: adv_ft(is_:ie_,js_:je_,nphi_),    var(is_:ie_,js_:je_,nphi_,3)
      real*8, intent(in) :: uvel(is_:ie_,js_:je_,nphi_),vvel(is_:ie_,js_:je_,nphi_),wvel(is_:ie_,js_:je_,nphi_)
      integer :: i,j,k,km1,kp2
      real*8 :: ddt,Rjp,Rj,Rjm,uCFL=0.5,Cr,Limiter
      Limiter(Cr)=max(0.D0,max(min(1.D0,2.D0*Cr), min(2.D0,Cr))) 

      do k=2,nphi-1
       do j=js_pe,je_pe
        do i=is_pe-1,ie_pe
         uCFL = ABS( uvel(i,j,k)*ddt/(cost(j)*dxt( min(nx,max(1,i)) )) )
         Rjp=(var(i+2,j,k,nn)-var(i+1,j,k,nn))*maskUp(i+1,j,k)
         Rj =(var(i+1,j,k,nn)-var(i  ,j,k,nn))*maskUp(i  ,j,k)
         Rjm=(var(i  ,j,k,nn)-var(i-1,j,k,nn))*maskUp(i-1,j,k)
         IF (Rj.NE.0.) THEN
          IF (uvel(i,j,k).GT.0) THEN; Cr=Rjm/Rj; ELSE; Cr=Rjp/Rj; ENDIF
         ELSE
          IF (uvel(i,j,k).GT.0) THEN; Cr=Rjm*1.E20; ELSE; Cr=Rjp*1.E20; ENDIF
         ENDIF
         Cr=Limiter(Cr)
         adv_fe(i,j,k) = uvel(i,j,k)*(var(i+1,j,k,nn)+var(i,j,k,nn))*0.5d0   &
                                -ABS(uvel(i,j,k))*((1.-Cr)+uCFL*Cr)*Rj*0.5d0
        enddo
       enddo
      enddo

      do k=2,nphi-1
       do j=js_pe-1,je_pe
        do i=is_pe,ie_pe
         Rjp=(var(i,j+2,k,nn)-var(i,j+1,k,nn))*maskVp(i,j+1,k)
         Rj =(var(i,j+1,k,nn)-var(i,j  ,k,nn))*maskVp(i,j  ,k)
         Rjm=(var(i,j  ,k,nn)-var(i,j-1,k,nn))*maskVp(i,j-1,k)
         uCFL = ABS( vvel(i,j,k)*ddt/dyt( min(ny,max(1,j)) ) )
         IF (Rj.NE.0.) THEN
          IF (vvel(i,j,k).GT.0) THEN; Cr=Rjm/Rj; ELSE; Cr=Rjp/Rj; ENDIF
         ELSE
          IF (vvel(i,j,k).GT.0) THEN; Cr=Rjm*1.E20; ELSE; Cr=Rjp*1.E20; ENDIF
         ENDIF
         Cr=Limiter(Cr)
         adv_fn(i,j,k) = vvel(i,j,k)*(var(i,j+1,k,nn)+var(i,j,k,nn))*0.5d0   &
                                -ABS(vvel(i,j,k))*((1.-Cr)+uCFL*Cr)*Rj*0.5d0
        enddo
       enddo
      enddo
 
      do k=1,nphi-1
       kp2=k+2; if (kp2>nphi) kp2=3
       km1=k-1; if (km1<1) km1=nphi-2
       do j=js_pe,je_pe
        do i=is_pe,ie_pe
         Rjp=(var(i,j,kp2,nn)-var(i,j,k+1,nn))*maskWp(i,j,k+1)
         Rj =(var(i,j,k+1,nn)-var(i,j,k  ,nn))*maskWp(i,j,k  )
         Rjm=(var(i,j,k  ,nn)-var(i,j,km1,nn))*maskWp(i,j,km1)
         uCFL = ABS( wvel(i,j,k)*ddt/dphit(k) )
         IF (Rj.NE.0.) THEN
          IF (wvel(i,j,k).GT.0) THEN; Cr=Rjm/Rj; ELSE; Cr=Rjp/Rj; ENDIF
         ELSE
          IF (wvel(i,j,k).GT.0) THEN; Cr=Rjm*1.E20; ELSE; Cr=Rjp*1.E20; ENDIF
         ENDIF
         Cr=Limiter(Cr)
         adv_ft(i,j,k) = wvel(i,j,k)*(var(i,j,k+1,nn)+var(i,j,k,nn))*0.5d0   &
                                -ABS(wvel(i,j,k))*((1.-Cr)+uCFL*Cr)*Rj*0.5d0
        enddo
       enddo
      enddo
end subroutine adv_flux_rossmix_superbee






subroutine adv_flux_rossmix_upwind(is_,ie_,js_,je_,nphi_,adv_fe,adv_fn,adv_ft,var,nn,dtloc,uvel,vvel,wvel)
!---------------------------------------------------------------------------------
! Calculates advection using upwind scheme
!---------------------------------------------------------------------------------
      use main_module   
      use rossmix_module   
      implicit none
      integer, intent(in) :: is_,ie_,js_,je_,nphi_,nn
      real*8, intent(inout) :: adv_fe(is_:ie_,js_:je_,nphi_), adv_fn(is_:ie_,js_:je_,nphi_)
      real*8, intent(inout) :: adv_ft(is_:ie_,js_:je_,nphi_),    var(is_:ie_,js_:je_,nphi_,3)
      real*8, intent(in) :: uvel(is_:ie_,js_:je_,nphi_),vvel(is_:ie_,js_:je_,nphi_),wvel(is_:ie_,js_:je_,nphi_)
      integer :: i,j,k
      real*8 :: Rj,dtloc

      do k=1,nphi
       do j=js_pe,je_pe
        do i=is_pe-1,ie_pe
         Rj =(var(i+1,j,k,nn)-var(i  ,j,k,nn))*maskUp(i  ,j,k)
         adv_fe(i,j,k) = uvel(i,j,k)*(var(i+1,j,k,nn)+var(i,j,k,nn))*0.5d0  -ABS(uvel(i,j,k))*Rj*0.5d0
        enddo
       enddo
      enddo

      do k=1,nphi
       do j=js_pe-1,je_pe
        do i=is_pe,ie_pe
         Rj =(var(i,j+1,k,nn)-var(i,j  ,k,nn))*maskVp(i,j  ,k)
         adv_fn(i,j,k) = cosu(j)*vvel(i,j,k)*(var(i,j+1,k,nn)+var(i,j,k,nn))*0.5d0 -ABS(cosu(j)*vvel(i,j,k))*Rj*0.5d0
        enddo
       enddo
      enddo
 
      do k=1,nphi-1
       do j=js_pe,je_pe
        do i=is_pe,ie_pe
         Rj =(var(i,j,k+1,nn)-var(i,j,k  ,nn))*maskWp(i,j,k  )
         adv_ft(i,j,k) = wvel(i,j,k)*(var(i,j,k+1,nn)+var(i,j,k,nn))*0.5d0 -ABS(wvel(i,j,k))*Rj*0.5d0
        enddo
       enddo
      enddo
      adv_ft(:,:,nz)=0.0
end subroutine adv_flux_rossmix_upwind








subroutine diff_flux_rossmix(is_,ie_,js_,je_,nphi_,fle,fln,flt,var,nn)
!---------------------------------------------------------------------------------
! Diffusion, taup1 is only updated but not set
!---------------------------------------------------------------------------------
 use main_module   
 use rossmix_module   
 implicit none
 integer :: is_,ie_,js_,je_,nphi_,nn
 real*8, dimension(is_:ie_,js_:je_,nphi_,3) :: var
 real*8, dimension(is_:ie_,js_:je_,nphi_) :: fle,fln,flt
 integer :: i,j,n
 do j=js_pe,je_pe
   do i=is_pe-1,ie_pe
    fle(i,j,:)=rossmix_lat_diffusion*(var(i+1,j,:,nn)-var(i,j,:,nn))/(cost(j)*dxu(i))*maskUp(i,j,:)
   enddo
 enddo
 do j=js_pe-1,je_pe
   do i=is_pe,ie_pe
    fln(i,j,:)=rossmix_lat_diffusion*(var(i,j+1,:,nn)-var(i,j,:,nn))/dyu(j)*maskVp(i,j,:)*cosu(j)
   enddo
 enddo
 do n=1,nphi-1
    flt(:,:,n)=rossmix_phi_diffusion*(var(:,:,n+1,nn)-var(:,:,n,nn))/dphiu(n)*maskWp(:,:,n)
 enddo
end subroutine diff_flux_rossmix




subroutine rossmix_check_cfl(m)
  use main_module   
  use rossmix_module   
  implicit none
  integer, intent(in) :: m
  integer :: i,j,n
  real*8 :: cfl,dt_loc
  logical :: stopit


  dt_loc = dt_tracer
  if (m==1) dt_loc = dt_tracer/rossmix_barotropic_fac

  stopit = .false.
  do n=2,nphi-1
   do j=js_pe,je_pe
     do i=is_pe,ie_pe
        cfl = abs( cgu_l(i,j,n,m)*dt_loc/(cost(j)*dxt(i)) )
        if (cfl>0.5) print*,' WARNING: CFL for cgu_l =',cfl,' at PE=',my_pe,' i=',i,' j=',j,' n=',n,' m=',m
        !if (cfl>0.5) exit
        if (cfl>0.5) stopit = .true.
         
        cfl = abs( cgu_s(i,j,n,m)*dt_loc/(cost(j)*dxt(i)) )
        if (cfl>0.5) print*,' WARNING: CFL for cgu_s =',cfl,' at PE=',my_pe,' i=',i,' j=',j,' n=',n,' m=',m
        if (cfl>0.5) stopit = .true.
        !if (cfl>0.5) exit
    enddo
    !if (cfl>0.5) exit
   enddo
   !if (cfl>0.5) exit
  enddo

  do n=2,nphi-1
   do j=js_pe,je_pe
    do i=is_pe,ie_pe
        cfl = abs( cgv_l(i,j,n,m)*dt_loc/(dyt(j))  )
        if (cfl>0.5) print*,' WARNING: CFL for cgv_l=',cfl,' at PE=',my_pe,' i=',i,' j=',j,' n=',n,' m=',m
        if (cfl>0.5) stopit = .true.
        !if (cfl>0.5) exit
        cfl = abs( cgv_s(i,j,n,m)*dt_loc/(dyt(j))  )
        if (cfl>0.5) print*,' WARNING: CFL for cgv_s=',cfl,' at PE=',my_pe,' i=',i,' j=',j,' n=',n,' m=',m
        if (cfl>0.5) stopit = .true.
        !if (cfl>0.5) exit
    enddo
     !if (cfl>0.5) exit
   enddo
     !if (cfl>0.5) exit
  enddo


  do n=2,nphi-1
   do j=js_pe,je_pe
    do i=is_pe,ie_pe
        cfl = abs( phidot_l(i,j,n,m)*dt_loc/dphit(n)  )
        if (cfl>0.5) print*,' WARNING: CFL for phidot_l =',cfl,' at PE=',my_pe,' i=',i,' j=',j,' n=',n,' m=',m
        if (cfl>0.5) stopit = .true.
        !if (cfl>0.5) exit
        cfl = abs( phidot_s(i,j,n,m)*dt_loc/dphit(n)  )
        if (cfl>0.5) print*,' WARNING: CFL for phidot_s =',cfl,' at PE=',my_pe,' i=',i,' j=',j,' n=',n,' m=',m
        if (cfl>0.5) stopit = .true.
        !if (cfl>0.5) exit
    enddo
     !if (cfl>0.5) exit
   enddo
     !if (cfl>0.5) exit
  enddo

  !if (stopit) call halt_stop('CFL violation')

end subroutine rossmix_check_cfl



