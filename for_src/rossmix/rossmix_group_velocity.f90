



 

subroutine rossmix_group_velocity(m)
  use main_module   
  use rossmix_module   
  implicit none
  integer, intent(in) :: m
  integer :: i,j,k,n
  real*8 :: Aloc(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx)
  real*8 :: Bloc(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx)
  real*8 :: Un(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx)
  real*8 :: Vn(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx)
  real*8 :: Ux(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx)
  real*8 :: Uy(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx)
  real*8 :: Vx(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx)
  real*8 :: Vy(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx)

  !   \int_{-h}^0  \v U phi_n^2 } dz 
  Un=0.0; Vn=0.0
  do k=1,nz
   do j=js_pe-onx,je_pe+onx-1
    do i=is_pe-onx,ie_pe+onx-1
       Un(i,j) = Un(i,j) + dzt(k)*0.5*( phin(i,j,k,m)**2 + phin(i+1,j,k,m)**2 )*u(i,j,k,tau)*maskU(i,j,k)
       Vn(i,j) = Vn(i,j) + dzt(k)*0.5*( phin(i,j,k,m)**2 + phin(i,j+1,k,m)**2 )*v(i,j,k,tau)*maskV(i,j,k)
    enddo
   enddo
  enddo

  call calc_gradients_ross(Un,Vn,Ux,Uy,Vx,Vy,is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx)

  if (enable_rossmix_lateral_stress) then 
   !---------------------------------------------------------------------------------
   ! growth rate nu
   !---------------------------------------------------------------------------------
   do n=2,nphi-1
     Aloc(:,:) = Ux(:,:)*cos(phit(n))**2 + Vy(:,:)*sin(phit(n))**2 + sin(phit(n))*cos(phit(n))*(Vx(:,:)+Uy(:,:))
     nu_l(:,:,n,m) = phi_nu_l(:,:,m)*Aloc(:,:)*maskTp(:,:,n)
     nu_s(:,:,n,m) = phi_nu_s(:,:,m)*Aloc(:,:)*maskTp(:,:,n)
   enddo
  endif

  do k=2,nphi-1
    cgu_l(:,:,k,m) = Un(:,:)*maskUp(:,:,k)
    cgv_l(:,:,k,m) = Vn(:,:)*maskVp(:,:,k)
    cgu_s(:,:,k,m) = Un(:,:)*maskUp(:,:,k)
    cgv_s(:,:,k,m) = Vn(:,:)*maskVp(:,:,k)
  enddo

 !  -  \beta e_1 cfac1 
  Aloc(:,:) = beta(:,:)*cfac1_l(:,:,m)
  Bloc(:,:) = beta(:,:)*cfac1_s(:,:,m)
  do i=is_pe-onx,ie_pe+onx-1
    Aloc(i,:) = (Aloc(i,:) + Aloc(i+1,:))*0.5*maskU(i,:,nz)
    Bloc(i,:) = (Bloc(i,:) + Bloc(i+1,:))*0.5*maskU(i,:,nz)
  enddo
  do k=2,nphi-1
    cgu_l(:,:,k,m) = cgu_l(:,:,k,m) - Aloc(:,:)*maskUp(:,:,k)
    cgu_s(:,:,k,m) = cgu_s(:,:,k,m) - Bloc(:,:)*maskUp(:,:,k)
  enddo

  ! + 2 \v n   n_1  beta  cfac2
  Aloc(:,:) = beta(:,:)*cfac2_l(:,:,m)
  Bloc(:,:) = beta(:,:)*cfac2_s(:,:,m)
  do i=is_pe-onx,ie_pe+onx-1
    Aloc(i,:) = (Aloc(i,:) + Aloc(i+1,:))*0.5*maskU(i,:,nz)
    Bloc(i,:) = (Bloc(i,:) + Bloc(i+1,:))*0.5*maskU(i,:,nz)
  enddo
  do k=2,nphi-1
     cgu_l(:,:,k,m) = cgu_l(:,:,k,m) + 2*cos(phit(k)) * cos(phit(k))*Aloc(:,:)*maskUp(:,:,k)
     cgu_s(:,:,k,m) = cgu_s(:,:,k,m) + 2*cos(phit(k)) * cos(phit(k))*Bloc(:,:)*maskUp(:,:,k)
  enddo

  Aloc(:,:) = beta(:,:)*cfac2_l(:,:,m)
  Bloc(:,:) = beta(:,:)*cfac2_s(:,:,m)
  do j=js_pe-onx,je_pe+onx-1
    Aloc(:,j) = (Aloc(:,j) + Aloc(:,j+1))*0.5*maskV(:,j,nz)
    Bloc(:,j) = (Bloc(:,j) + Bloc(:,j+1))*0.5*maskV(:,j,nz)
  enddo
  do k=2,nphi-1
     cgv_l(:,:,k,m) = cgv_l(:,:,k,m) + 2*sin(phit(k)) * cos(phit(k))*Aloc(:,:)*maskVp(:,:,k)
     cgv_s(:,:,k,m) = cgv_s(:,:,k,m) + 2*sin(phit(k)) * cos(phit(k))*Bloc(:,:)*maskVp(:,:,k)
  enddo


  !   \v n \cdot \rvec{\vn} ( \v U_n \cdot \v n )  
  ! simpler with
  !  \v n \cdot \rvec{\vn} (\v U \cdot \v n) = 
  !   \sin^2 \phi \px V  - \cos^2 \phi \py U   + \sin \phi \cos \phi (\px U  -  \py V) 

  phidot_s(:,:,:,m) = 0.; phidot_l(:,:,:,m)= 0.
  do n=1,nphi
   Aloc(:,:) = sin(phit(n))**2*Vx(:,:)-cos(phit(n))**2*Uy(:,:)+sin(phit(n))*cos(phit(n))*(Ux(:,:)-Vy(:,:))
   phidot_l(:,:,n,m) = phidot_l(:,:,n,m) + Aloc(:,:)
   phidot_s(:,:,n,m) = phidot_s(:,:,n,m) + Aloc(:,:)
  enddo
  phidot_l(:,:,:,m) = phidot_l(:,:,:,m)*maskWp
  phidot_s(:,:,:,m) = phidot_s(:,:,:,m)*maskWp

  call border_exchg_xyp(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nphi,cgu_l(:,:,:,m)); 
  call setcyclic_xyp(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nphi,cgu_l(:,:,:,m))
  call border_exchg_xyp(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nphi,cgu_s(:,:,:,m)); 
  call setcyclic_xyp(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nphi,cgu_s(:,:,:,m))
  call border_exchg_xyp(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nphi,cgv_l(:,:,:,m)); 
  call setcyclic_xyp(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nphi,cgv_l(:,:,:,m))
  call border_exchg_xyp(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nphi,cgv_s(:,:,:,m)); 
  call setcyclic_xyp(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nphi,cgv_s(:,:,:,m))
  call border_exchg_xyp(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nphi,phidot_l(:,:,:,m)); 
  call setcyclic_xyp(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nphi,phidot_l(:,:,:,m))
  call border_exchg_xyp(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nphi,phidot_s(:,:,:,m)); 
  call setcyclic_xyp(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nphi,phidot_s(:,:,:,m))

end subroutine rossmix_group_velocity





subroutine rossmix_growth_rate(m)
 use main_module   
 use rossmix_module   
 implicit none
 integer :: n,m,k,j,i
 real*8 :: Aloc(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx),fxa,fxb,fxc
 real*8 :: Bloc(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx)

 !---------------------------------------------------------------------------------
 ! growth rate gamma_n : long wave compartment
 !---------------------------------------------------------------------------------
  Aloc = 0.; Bloc = 0.
  do k=1,nz-1
     do j=js_pe-onx,je_pe+onx-1
      do i=is_pe-onx,ie_pe+onx-1
       fxa = 0.5*(phi_ga_l(i,j,k,m) +phi_ga_l(i+1,j,k,m)) 
       fxb = sqrt( max(rossmix_N2min, 0.5*(Nsqr(i,j,k,tau) +Nsqr(i+1,j,k,tau))  ) )
       fxc =  (u(i,j,k+1,tau) - u(i,j,k,tau))/dzw(k)
       Aloc(i,j) = Aloc(i,j) + dzw(k)*maskU(i,j,k)*maskU(i,j,k+1)*fxa*abs(coriolis_t(i,j))/fxb*fxc
       fxa = 0.5*(phi_ga_l(i,j,k,m) +phi_ga_l(i,j+1,k,m))
       fxb = sqrt( max(rossmix_N2min, 0.5*(Nsqr(i,j,k,tau) +Nsqr(i,j+1,k,tau)) ))
       fxc =  (v(i,j,k+1,tau) - v(i,j,k,tau))/dzw(k)
       Bloc(i,j) = Bloc(i,j) + dzw(k)*maskV(i,j,k)*maskV(i,j,k+1)*fxa*abs(coriolis_t(i,j))/fxb*fxc
      enddo
     enddo
  enddo
  ! interpolate on T grid 
  call ugrid_to_tgrid( is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,1,Aloc,Aloc)  
  call vgrid_to_tgrid( is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,1,Bloc,Bloc)  
  do n=1,nphi
      gamma_l(:,:,n,m) = (cos(phit(n))*Aloc(:,:) + sin(phit(n))*Bloc(:,:))/max(1D-12,Rn(:,:,m))
  enddo
 !---------------------------------------------------------------------------------
 ! growth rate gamma_n : short wave compartment
 !---------------------------------------------------------------------------------
  Aloc = 0.; Bloc = 0.
  do k=1,nz-1
     do j=js_pe-onx,je_pe+onx-1
      do i=is_pe-onx,ie_pe+onx-1
       fxa = 0.5*(phi_ga_s(i,j,k,m) +phi_ga_s(i+1,j,k,m)) 
       fxb = sqrt( max(rossmix_N2min, 0.5*(Nsqr(i,j,k,tau) +Nsqr(i+1,j,k,tau)) ))
       fxc =  (u(i,j,k+1,tau) - u(i,j,k,tau))/dzw(k)
       Aloc(i,j) = Aloc(i,j) + dzw(k)*maskU(i,j,k)*maskU(i,j,k+1)*fxa*abs(coriolis_t(i,j))/fxb*fxc
       fxa = 0.5*(phi_ga_s(i,j,k,m) +phi_ga_s(i,j+1,k,m))
       fxb = sqrt( max(rossmix_N2min, 0.5*(Nsqr(i,j,k,tau) +Nsqr(i,j+1,k,tau)) ))
       fxc =  (v(i,j,k+1,tau) - v(i,j,k,tau))/dzw(k)
       Bloc(i,j) = Bloc(i,j) + dzw(k)*maskV(i,j,k)*maskV(i,j,k+1)*fxa*abs(coriolis_t(i,j))/fxb*fxc
      enddo
     enddo
  enddo
  ! interpolate on T grid 
  call ugrid_to_tgrid( is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,1,Aloc,Aloc)  
  call vgrid_to_tgrid( is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,1,Bloc,Bloc)  
  do n=1,nphi
      gamma_s(:,:,n,m) = (cos(phit(n))*Aloc(:,:) + sin(phit(n))*Bloc(:,:))/max(1D-12,Rn(:,:,m))
  enddo


 if (rossmix_symmetrisation > 0d0 ) then
   fxa = rossmix_symmetrisation
   lambda_l(:,:,:,m) = rossmix_c_lambda*(sqrt(4*gamma_l(:,:,:,m)**2+fxa**2/4)-fxa/2.)
   lambda_s(:,:,:,m) = rossmix_c_lambda*(sqrt(4*gamma_s(:,:,:,m)**2+fxa**2/4)-fxa/2.)

 elseif (rossmix_c_tau > 0d0 ) then
   do i=is_pe-onx,ie_pe+onx
    do j=js_pe-onx,je_pe+onx
      fxa = rossmix_c_tau*beta(i,j)*Rn(i,j,m)
      lambda_l(i,j,:,m) = rossmix_c_lambda*(sqrt(4*gamma_l(i,j,:,m)**2+fxa**2/4)-fxa/2.)
      lambda_s(i,j,:,m) = rossmix_c_lambda*(sqrt(4*gamma_s(i,j,:,m)**2+fxa**2/4)-fxa/2.)
    enddo  
   enddo  
 elseif (enable_rossmix_use_gamma) then
   lambda_l(:,:,:,m) = 2*gamma_l(:,:,:,m)
   lambda_s(:,:,:,m) = 2*gamma_s(:,:,:,m)
 endif

end subroutine rossmix_growth_rate







subroutine calc_gradients_ross(Uloc,Vloc,Uxloc,Uyloc,Vxloc,Vyloc,is_,ie_,js_,je_)

! return Ux,Uy,Vx,Vy on T grid

  use main_module   
  implicit none
  integer :: is_,ie_,js_,je_
  real*8, dimension(is_:ie_,js_:je_) :: Uloc,Vloc,Uxloc,Uyloc,Vxloc,Vyloc
  integer :: i,j
  real*8 :: aloc(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx),fxa

  do j=js_pe-onx,je_pe+onx
   do i=is_pe-onx+1,ie_pe+onx
     Uxloc(i,j) = (Uloc(i,j) - Uloc(i-1,j))/(cost(j)*dxt(i)) *maskT(i,j,nz)
   enddo
  enddo

  do j=js_pe-onx+1,je_pe+onx
    Vyloc(:,j) = (Vloc(:,j) - Vloc(:,j-1))/dyt(j)*maskT(:,j,nz)
  enddo

  do j=js_pe-onx,je_pe+onx-1
   Aloc(:,j) = ( Uloc(:,j+1) - Uloc(:,j) )/dyu(j) *maskZ(:,j,nz)
  enddo
  Aloc(:,je_pe+onx) =  Aloc(:,je_pe+onx-1)
  do j=js_pe-onx+1,je_pe+onx
   do i=is_pe-onx+1,ie_pe+onx
     fxa = max(1d-12,maskZ(i,j,nz)+maskZ(i-1,j,nz)+maskZ(i,j-1,nz)+maskZ(i-1,j-1,nz) )
     Uyloc(i,j) = (Aloc(i,j)+Aloc(i-1,j)+Aloc(i,j-1)+Aloc(i-1,j-1))/fxa*maskT(i,j,nz)
   enddo
  enddo

  do j=js_pe-onx,je_pe+onx
   do i=is_pe-onx,ie_pe+onx-1
    Aloc(i,j) = ( Vloc(i+1,j) - Vloc(i,j) )/(cosu(j)*dxu(i))*maskZ(i,j,nz)
   enddo
  enddo
  Aloc(ie_pe+onx,:) =  Aloc(ie_pe+onx-1,:)
  do j=js_pe-onx+1,je_pe+onx
   do i=is_pe-onx+1,ie_pe+onx
     fxa = max(1d-12,maskZ(i,j,nz)+maskZ(i-1,j,nz)+maskZ(i,j-1,nz)+maskZ(i-1,j-1,nz))
     Vxloc(i,j) = (Aloc(i,j)+Aloc(i-1,j)+Aloc(i,j-1)+Aloc(i-1,j-1))/fxa *maskT(i,j,nz)
   enddo
  enddo


end subroutine calc_gradients_ross


