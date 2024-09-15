


subroutine rossmix2_main
 use main_module   
 use rossmix2_module   
 use diagnostics_module 
 implicit none
 logical, save  :: first = .true.
 integer :: i,j,k,n
 real*8 :: fxa,get_rho
 
 if (.not. enable_rossmix2) return

 if (.not. enable_thermodynamic_equation ) then
 
  do k=1,nz
   do j=js_pe-onx,je_pe+onx
    do i=is_pe-onx,ie_pe+onx
     rho(i,j,k,tau) = get_rho(salt(i,j,k,tau),temp(i,j,k,tau),abs(zt(k)))*maskT(i,j,k)
    enddo
   enddo
  enddo
 
  do k=1,nz-1
   do j=js_pe-onx,je_pe+onx
    do i=is_pe-onx,ie_pe+onx
     fxa =  -grav/rho_0/dzw(k)*maskW(i,j,k)
     Nsqr(i,j,k,tau) =fxa*(get_rho(salt(i,j,k+1,tau),temp(i,j,k+1,tau),abs(zt(k)))-rho(i,j,k,tau)) 
    enddo
   enddo
  enddo
  Nsqr(:,:,nz,tau)=Nsqr(:,:,nz-1,tau)
 endif

 ! limit Nsqr for later use
 Nsqr_lim =  max(rossmix2_N2min,Nsqr(:,:,:,tau))

 if (first .or. mod(itt,int(rossmix2_calc_parameter_period/dt_tracer))  == 0) then
  call rossmix2_vertical_mode
  call rossmix2_group_velocity
  if (enable_rossmix2_simple) then
   call rossmix2_eigenvalue_simple
  else if (enable_rossmix2_lapack) then
   call rossmix2_eigenvalue_lapack
  else 
   call rossmix2_eigenvalue_approx
  endif
  !
 endif
 
 if (itt==0 .or. mod(itt,int(snapint/dt_tracer))  == 0)  call rossmix2_diag_snap
 
  
 ! B = \int dphi \rvec{u} b /N^2 = - \int dphi \v n b /N^2
 Bx=0
 By=0
 do n=2,nphi-1
  do k=1,nz-1
   do j=js_pe,je_pe
    do i=is_pe,ie_pe
     fxa = ub(i,j,k,n)/Nsqr_lim(i,j,k)*maskW(i,j,k) 
     Bx(i,j,k) = Bx(i,j,k) - E_r(i,j,n,tau)*fxa*dphit(n)*coriolis_t(i,j)*cos(phit(n))
     By(i,j,k) = By(i,j,k) - E_r(i,j,n,tau)*fxa*dphit(n)*coriolis_t(i,j)*sin(phit(n))     
    enddo
   enddo
  enddo
 enddo 

 call border_exchg_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,Bx) 
 call setcyclic_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,Bx)
 call border_exchg_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,By) 
 call setcyclic_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,By)

 call tgrid_to_ugrid( is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,By,By)  
 call tgrid_to_vgrid( is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,Bx,Bx)  
  
 k=1
 do j=js_pe,je_pe
  do i=is_pe,ie_pe
     ue(i,j,k) = -(By(i,j,k))/dzt(k)*maskU(i,j,k)
     ve(i,j,k) =  (Bx(i,j,k))/dzt(k)*maskV(i,j,k)
  enddo
 enddo
 do k=2,nz
  do j=js_pe,je_pe
    do i=is_pe,ie_pe
     ue(i,j,k) = -(By(i,j,k)-By(i,j,k-1))/dzt(k)*maskU(i,j,k)
     ve(i,j,k) =  (Bx(i,j,k)-Bx(i,j,k-1))/dzt(k)*maskV(i,j,k)
    enddo
  enddo
 enddo

 call border_exchg_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,ue) 
 call setcyclic_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,ue)
 call border_exchg_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,ve) 
 call setcyclic_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,ve)

 k=1
 do j=js_pe-onx+1,je_pe+onx
   do i=is_pe-onx+1,ie_pe+onx
     we(i,j,k) =-maskW(i,j,k)*dzt(k)* &
               ((        ue(i,j,k)-          ue(i-1,j,k))/(cost(j)*dxt(i)) &
               +(cosu(j)*ve(i,j,k)-cosu(j-1)*ve(i,j-1,k))/(cost(j)*dyt(j)) )
   enddo
 enddo
 do k=2,nz
  do j=js_pe-onx+1,je_pe+onx
   do i=is_pe-onx+1,ie_pe+onx
     we(i,j,k) = we(i,j,k-1)-maskW(i,j,k)*dzt(k)* &
               ((        ue(i,j,k)          -ue(i-1,j,k))/(cost(j)*dxt(i)) &
               +(cosu(j)*ve(i,j,k)-cosu(j-1)*ve(i,j-1,k))/(cost(j)*dyt(j)) )
   enddo
  enddo
 enddo 
  
  
 ! lateral stresses
 stress_ux = 0; stress_uy = 0; stress_vx = 0; stress_vy = 0 
 do n=2,nphi-1  
   do k=1,nz
    stress_ux(:,:,k) = stress_ux(:,:,k) + dphit(n)*sin(phit(n))*sin(phit(n))*uu(:,:,k,n)*E_r(:,:,n,tau)
    stress_uy(:,:,k) = stress_uy(:,:,k) - dphit(n)*sin(phit(n))*cos(phit(n))*uu(:,:,k,n)*E_r(:,:,n,tau) 
    stress_vx(:,:,k) = stress_vx(:,:,k) - dphit(n)*sin(phit(n))*cos(phit(n))*uu(:,:,k,n)*E_r(:,:,n,tau) 
    stress_vy(:,:,k) = stress_vy(:,:,k) + dphit(n)*cos(phit(n))*cos(phit(n))*uu(:,:,k,n)*E_r(:,:,n,tau)
   enddo
 enddo

 !---------------------------------------------------------------------------------
 ! Interpolation
 !---------------------------------------------------------------------------------
 stress_ux(is_pe-onx:ie_pe+onx-1,:,:) = stress_ux(is_pe-onx+1:ie_pe+onx,:,:)
 call tgrid_to_ugrid( is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,stress_uy,stress_uy)  
 call tgrid_to_vgrid( is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,stress_uy,stress_uy)  

 stress_vy(is_pe-onx:ie_pe+onx-1,:,:) = stress_vy(is_pe-onx+1:ie_pe+onx,:,:)
 call tgrid_to_ugrid( is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,stress_vx,stress_vx)  
 call tgrid_to_vgrid( is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,stress_vx,stress_vx)  

 !---------------------------------------------------------------------------------
 ! Masking
 !---------------------------------------------------------------------------------
 do i=is_pe-onx,ie_pe+onx-1
   stress_ux(i,:,:) = stress_ux(i,:,:)*maskU(i,:,:)*maskU(i+1,:,:)
 enddo
 do j=js_pe-onx,je_pe+onx-1
    stress_uy(:,j,:) = stress_uy(:,j,:)*maskU(:,j,:)*maskU(:,j+1,:)
 enddo
 do i=is_pe-onx,ie_pe+onx-1
   stress_vx(i,:,:) = stress_vx(i,:,:)*maskV(i,:,:)*maskV(i+1,:,:)
 enddo
 do j=js_pe-onx,je_pe+onx-1
   stress_vy(:,j,:) = stress_vy(:,j,:)*maskV(:,j,:)*maskV(:,j+1,:)
 enddo
   
 first = .false.
end subroutine rossmix2_main



subroutine rossmix2_initialize
  use main_module   
  use rossmix2_module   
  implicit none
  integer :: i,j,k

  if (.not. enable_rossmix2) return
 
  if (enable_rossmix2_include_lateral_friction) then
    if (.not.enable_conserve_energy) then
     if (my_pe==0) print*,' conflict with lateral dissipation redirection in Rossmix2 '
     call halt_stop(' in rossmix2_initialize (3)')
    endif    
    if (.not.enable_store_lateral_friction_heat) then
     if (my_pe==0) print*,' conflict of lateral dissipation redirection in Rossmix2 '
     call halt_stop(' in rossmix2_initialize (4)')
    endif
  endif
  
  ! wavenumber grid  
  dphit=2.*pi/(nphi-2);  dphiu=dphit
  phit(1)=0.0-dphit(1)/2.; phiu(1)=phit(1)+dphit(1)/2.
  do i=2,nphi
   phit(i)=phit(i-1)+dphit(i); phiu(i)=phiu(i-1)+dphiu(i)
  enddo

  ! topographic mask for waves
  maskTp=0.0
  do j=js_pe,je_pe
    do i=is_pe,ie_pe
      if ( kbot(i,j) /=0 ) maskTp(i,j,:)=1.0
    enddo
  enddo
  call border_exchg_xyp(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nphi,maskTp) 
  call setcyclic_xyp   (is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nphi,maskTp)
  maskUp=maskTp
  do i=is_pe-onx,ie_pe+onx-1
    maskUp(i,:,:)=min(maskTp(i,:,:),maskTp(i+1,:,:))
  enddo
  call border_exchg_xyp(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nphi,maskUp)
  call setcyclic_xyp   (is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nphi,maskUp)
  maskVp=maskTp
  do j=js_pe-onx,je_pe+onx-1
    maskVp(:,j,:)=min(maskTp(:,j,:),maskTp(:,j+1,:))
  enddo
  call border_exchg_xyp(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nphi,maskVp)
  call setcyclic_xyp   (is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nphi,maskVp)
  maskWp=maskTp
  do k=1,nphi-1
     maskWp(:,:,k)=min(maskTp(:,:,k),maskTp(:,:,k+1))
  enddo
  
  call rossmix2_init_snap_cdf
end subroutine rossmix2_initialize




subroutine rossmix2_integrate
 use main_module   
 use rossmix2_module   
 implicit none
 integer :: i,j,n,k,np2,nm1
 real*8 :: aloc(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx)
 real*8 :: advfe(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nphi) 
 real*8 :: advfn(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nphi) 
 real*8 :: advfp(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nphi)
 real*8, dimension(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx) :: Ux,Uy,Vx,Vy

 real*8 :: Rjp,Rj,Rjm,uCFL=0.5,Cr,Limiter,fxa
 Limiter(Cr)=max(0.D0,max(min(1.D0,2.D0*Cr), min(2.D0,Cr))) 

 if (.not. enable_rossmix2) return

 ! advective and diffusive zonal flux of E 
 do n=2,nphi-1
  do j=js_pe,je_pe
   do i=is_pe-1,ie_pe
    uCFL = ABS( cgu(i,j,n)*dt_tracer/(cost(j)*dxt( min(nx,max(1,i)) )) )
    Rjp=(E_r(i+2,j,n,tau)-E_r(i+1,j,n,tau))*maskUp(i+1,j,n)
    Rj =(E_r(i+1,j,n,tau)-E_r(i  ,j,n,tau))*maskUp(i  ,j,n)
    Rjm=(E_r(i  ,j,n,tau)-E_r(i-1,j,n,tau))*maskUp(i-1,j,n)
    IF (Rj.NE.0.) THEN
          IF (cgu(i,j,n).GT.0) THEN; Cr=Rjm/Rj; ELSE; Cr=Rjp/Rj; ENDIF
    ELSE
          IF (cgu(i,j,n).GT.0) THEN; Cr=Rjm*1.E20; ELSE; Cr=Rjp*1.E20; ENDIF
    ENDIF
    Cr=Limiter(Cr)
    advfe(i,j,n) = -(cgu(i,j,n)*(E_r(i+1,j,n,tau)+E_r(i,j,n,tau))*0.5d0   &
                    -ABS(cgu(i,j,n))*((1.-Cr)+uCFL*Cr)*Rj*0.5d0 )* maskUp(i,j,n) &                 
          + rossmix2_Kh*(E_r(i+1,j,n,tau)-E_r(i,j,n,tau))/(cost(j)*dxu(i))*maskUp(i,j,n)              
   enddo
  enddo
 enddo

 ! advective and diffusive meridional flux of E
 do n=2,nphi-1
  do j=js_pe-1,je_pe
   do i=is_pe,ie_pe
    Rjp=(E_r(i,j+2,n,tau)-E_r(i,j+1,n,tau))*maskVp(i,j+1,n)
    Rj =(E_r(i,j+1,n,tau)-E_r(i,j  ,n,tau))*maskVp(i,j  ,n)
    Rjm=(E_r(i,j  ,n,tau)-E_r(i,j-1,n,tau))*maskVp(i,j-1,n)
    uCFL = ABS( cgv(i,j,n)*dt_tracer/dyt( min(ny,max(1,j)) ) )
    IF (Rj.NE.0.) THEN
          IF (cgv(i,j,n).GT.0) THEN; Cr=Rjm/Rj; ELSE; Cr=Rjp/Rj; ENDIF
    ELSE
          IF (cgv(i,j,n).GT.0) THEN; Cr=Rjm*1.E20; ELSE; Cr=Rjp*1.E20; ENDIF
    ENDIF
    Cr=Limiter(Cr)
    advfn(i,j,n) = -(cosu(j)*cgv(i,j,n)*(E_r(i,j+1,n,tau)+E_r(i,j,n,tau))*0.5d0   &
                -ABS(cosu(j)*cgv(i,j,n))*((1.-Cr)+uCFL*Cr)*Rj*0.5d0 ) * maskVp(i,j,n) &             
             + rossmix2_Kh*(E_r(i,j+1,n,tau)-E_r(i,j,n,tau))/dyu(j)*maskVp(i,j,n)*cosu(j)          
   enddo
  enddo
 enddo
 
! advective and diffusive flux in wave angle of E
 if (enable_rossmix2_calc_phidot) then
  do n=1,nphi-1
   np2=k+2; if (np2>nphi) np2=3
   nm1=n-1; if (nm1<1) nm1=nphi-2
   do j=js_pe,je_pe
    do i=is_pe,ie_pe
     Rjp=(E_r(i,j,np2,tau)-E_r(i,j,n+1,tau))*maskWp(i,j,n+1)
     Rj =(E_r(i,j,n+1,tau)-E_r(i,j,n  ,tau))*maskWp(i,j,n  )
     Rjm=(E_r(i,j,n  ,tau)-E_r(i,j,nm1,tau))*maskWp(i,j,nm1)
     uCFL = ABS( phidot(i,j,n)*dt_tracer/dphit(n) )
     IF (Rj.NE.0.) THEN
          IF (phidot(i,j,n).GT.0) THEN; Cr=Rjm/Rj; ELSE; Cr=Rjp/Rj; ENDIF
     ELSE
          IF (phidot(i,j,n).GT.0) THEN; Cr=Rjm*1.E20; ELSE; Cr=Rjp*1.E20; ENDIF
     ENDIF
     Cr=Limiter(Cr)
     advfp(i,j,n) = phidot(i,j,n)*(E_r(i,j,n+1,tau)+E_r(i,j,n,tau))*0.5d0   &
                              -ABS(phidot(i,j,n))*((1.-Cr)+uCFL*Cr)*Rj*0.5d0 &
                      + rossmix2_Kp*(E_r(i,j,n+1,tau)-E_r(i,j,n,tau))/dphiu(n)                         
    enddo
   enddo
  enddo
 else
  do n=1,nphi-1
   do j=js_pe,je_pe
    do i=is_pe,ie_pe
     advfp(i,j,n) =  rossmix2_Kp*(E_r(i,j,n+1,tau)-E_r(i,j,n,tau))/dphiu(n)                         
    enddo
   enddo
  enddo
 endif

 ! energy transfers forcing is +n_1 U_z + n_2 V_z
 forcx=0;forcy=0
 do n=2,nphi-1
  do k=1,nz-1  
   do j=js_pe-onx,je_pe+onx
    do i=is_pe-onx,ie_pe+onx
     fxa = ub(i,j,k,n)/Nsqr_lim(i,j,k) 
     aloc(i,j) = E_r(i,j,n,tau)*fxa*coriolis_t(i,j)**2*maskW(i,j,k) 
    enddo
   enddo
   call tgrid_to_ugrid_at_k(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,k,aloc,aloc)  
   do i=is_pe-1,ie_pe
    forcx(i,:,n) = forcx(i,:,n) + aloc(i,:)*cos(phit(n))*(u(i,:,k+1,tau)-u(i,:,k,tau)) &
                                 *maskW(i,:,k)*maskW(i+1,:,k)
   enddo
  enddo
  
  do k=1,nz-1  
   do j=js_pe-onx,je_pe+onx
    do i=is_pe-onx,ie_pe+onx
     fxa = ub(i,j,k,n)/Nsqr_lim(i,j,k)
     aloc(i,j) = E_r(i,j,n,tau)*fxa*coriolis_t(i,j)**2*maskW(i,j,k) 
    enddo
   enddo
   call tgrid_to_vgrid_at_k( is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,k,aloc,aloc)  
   do j=js_pe-1,je_pe
    forcy(:,j,n) = forcy(:,j,n) + aloc(:,j)*sin(phit(n))*(v(:,j,k+1,tau)-v(:,j,k,tau)) &
                                 *maskW(:,j,k)*maskW(:,j+1,k)
   enddo
  enddo

  ! interpolate ugrid to tgrid
  call ugrid_to_tgrid( is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,1,forcx(:,:,n),forcx(:,:,n))  
  call vgrid_to_tgrid( is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,1,forcy(:,:,n),forcy(:,:,n))  

 enddo

 ! lateral stress forcing 
 forcuu=0d0
 do n=2,nphi-1
  do k=1,nz
    call calc_gradients_ross(u(:,:,k,tau),v(:,:,k,tau),Ux,Uy,Vx,Vy,is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx) 
    forcuu(:,:,n) = forcuu(:,:,n) + maskTp(:,:,n)*dzt(k)*E_r(:,:,n,tau)*uu(:,:,k,n) & 
                                      *( Ux(:,:)*cos(phit(n))**2 + Vy(:,:)*sin(phit(n))**2 &
                                      + sin(phit(n))*cos(phit(n))*(Vx(:,:)+Uy(:,:))  ) 
  enddo
 enddo
  
  
 forc_diss = 0d0
 if (enable_rossmix2_include_lateral_friction) then
  do k=1,nz
    forc_diss(:,:) = forc_diss(:,:) + dzt(k)*K_diss_h(:,:,k)*maskT(:,:,k)*rossmix2_lateral_friction_fraction
  enddo
 endif

 
 ! E_n = E_(n-1) + dt*fxa - dt*E_n*E_(n-1)*a
 ! E_n  = (E_(n-1) + dt*fxa )/(1+ dt*E_(n-1)*a)

 do n=2,nphi-1   
  if (enable_rossmix2_linear_damping) then
   do j=js_pe,je_pe
    do i=is_pe,ie_pe
     fxa = maskTp(i,j,n)*( (advfe(i,j,n)- advfe(i-1,j,n))/(cost(j)*dxt(i)) &
                          +(advfn(i,j,n)- advfn(i,j-1,n))/(cost(j)*dyt(j)) &
                          +(advfp(i,j,n)- advfp(i,j,n-1))/dphit(n) )  
     fxa =  fxa + forc_diss(i,j)*dphit(n) + forcx(i,j,n) + forcy(i,j,n) + forcuu(i,j,n)                                        
     E_r(i,j,n,taup1) = E_r(i,j,n,tau) + dt_tracer*fxa - dt_tracer*E_r(i,j,n,tau)*rossmix2_damping 
    enddo
   enddo  
  else ! quadratic damping
   do j=js_pe,je_pe
    do i=is_pe,ie_pe
     fxa = maskTp(i,j,n)*( (advfe(i,j,n)- advfe(i-1,j,n))/(cost(j)*dxt(i)) &
                          +(advfn(i,j,n)- advfn(i,j-1,n))/(cost(j)*dyt(j)) &
                          +(advfp(i,j,n)- advfp(i,j,n-1))/dphit(n) )  
     fxa =  fxa + forc_diss(i,j)*dphit(n) + forcx(i,j,n) + forcy(i,j,n) + forcuu(i,j,n)
     E_r(i,j,n,taup1) = (E_r(i,j,n,tau) + dt_tracer*fxa)/(1+dt_tracer*E_r(i,j,n,tau)*rossmix2_damping)
    enddo
   enddo
  endif 
  
  if (enable_rossmix2_limit_energy)  E_r(:,:,n,taup1) = max(1d-12,E_r(:,:,n,taup1) )
  
 enddo


 !---------------------------------------------------------------------------------
 ! boundary exchange over PEs and cyclic boundary condition
 !---------------------------------------------------------------------------------
 call border_exchg_xyp(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nphi,E_r(:,:,:,taup1)) 
 call setcyclic_xyp(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nphi,E_r(:,:,:,taup1))

end subroutine rossmix2_integrate







