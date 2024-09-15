

subroutine rossmix2_group_velocity
 ! calculate group velocity, etc. for eddy fluxes
 use main_module   
 use rossmix2_module   
 implicit none
 integer :: i,j,k
 real*8, dimension(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx) :: Aloc,Ux,Uy,Vx,Vy

 !   \int_{-h}^0  \v U phi_1^2 } dz 
 U1=0.0; V1=0.0
 do k=1,nz
   do j=js_pe-onx,je_pe+onx-1
    do i=is_pe-onx,ie_pe+onx-1
       U1(i,j) = U1(i,j) + dzt(k)*0.5*( phi1(i,j,k)**2 + phi1(i+1,j,k)**2 )*u(i,j,k,tau)*maskU(i,j,k)
       V1(i,j) = V1(i,j) + dzt(k)*0.5*( phi1(i,j,k)**2 + phi1(i,j+1,k)**2 )*v(i,j,k,tau)*maskV(i,j,k)
    enddo
   enddo
 enddo

 do k=2,nphi-1
  cgu(:,:,k) = U1(:,:)*maskUp(:,:,k)
  cgv(:,:,k) = V1(:,:)*maskVp(:,:,k)
 enddo 
  
 !  -  \beta e_1 A4
 Aloc(:,:) = beta(:,:)*rossmix2_A4*R1(:,:)**2
 do i=is_pe-onx,ie_pe+onx-1
   Aloc(i,:) = (Aloc(i,:) + Aloc(i+1,:))*0.5*maskU(i,:,nz)
 enddo
 do k=2,nphi-1 
   cgu(:,:,k) = cgu(:,:,k) - Aloc(:,:)*maskUp(:,:,k)
 enddo

 ! + 2 \v n   n_1  beta A3
 Aloc(:,:) = beta(:,:)*rossmix2_A3*R1(:,:)**2
 do i=is_pe-onx,ie_pe+onx-1
    Aloc(i,:) = (Aloc(i,:) + Aloc(i+1,:))*0.5*maskU(i,:,nz)
 enddo
 do k=2,nphi-1
     cgu(:,:,k) = cgu(:,:,k) + 2*cos(phit(k)) * cos(phit(k))*Aloc(:,:)*maskUp(:,:,k)
 enddo

 Aloc(:,:) = beta(:,:)*rossmix2_A3*R1(:,:)**2
 do j=js_pe-onx,je_pe+onx-1
    Aloc(:,j) = (Aloc(:,j) + Aloc(:,j+1))*0.5*maskV(:,j,nz)
 enddo
 do k=2,nphi-1
     cgv(:,:,k) = cgv(:,:,k) + 2*sin(phit(k)) * cos(phit(k))*Aloc(:,:)*maskVp(:,:,k)
 enddo

 call border_exchg_xyp(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nphi,cgu(:,:,:)); 
 call setcyclic_xyp(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nphi,cgu(:,:,:))
 call border_exchg_xyp(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nphi,cgv(:,:,:)); 
 call setcyclic_xyp(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nphi,cgv(:,:,:))


 if (enable_rossmix2_calc_phidot) then
 
  !   \v n \cdot \rvec{\vn} ( \v U_n \cdot \v n )  
  ! simpler with
  !  \v n \cdot \rvec{\vn} (\v U \cdot \v n) = 
  !   \sin^2 \phi \px V  - \cos^2 \phi \py U   + \sin \phi \cos \phi (\px U  -  \py V) 
  call calc_gradients_ross(U1,V1,Ux,Uy,Vx,Vy,is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx)
  phidot = 0.
  do k=1,nphi
   phidot(:,:,k) = (sin(phit(k))**2*Vx(:,:) - cos(phit(k))**2*Uy(:,:) + &
                    sin(phit(k))*cos(phit(k))*(Ux(:,:)-Vy(:,:)) )*maskWp(:,:,k)
  enddo
  call border_exchg_xyp(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nphi,phidot)
  call setcyclic_xyp(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nphi,phidot)
 endif

end subroutine rossmix2_group_velocity






subroutine rossmix2_vertical_mode
  ! calculate approximate 1. vertical mode
  ! Nsqr and frequencies must be set before this routine
  use main_module   
  use rossmix2_module   
  implicit none
  integer :: k
  real*8 :: norm(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx)
  real*8 :: xi(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz)
  real*8 :: xiW(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz)
  real*8 :: c1(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx) 
  real*8 :: c1W(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx) 
  real*8 :: Nsqrt(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz)
  
  ! interpolate Nsqr 
  Nsqrt(:,:,1) = Nsqr_lim(:,:,1)
  do k=2,nz
    Nsqrt(:,:,k) = (Nsqr_lim(:,:,k)*maskW(:,:,k)+Nsqr_lim(:,:,k-1)*maskW(:,:,k-1)) / &
                                (maskW(:,:,k)+maskW(:,:,k-1)+1d-22)*maskT(:,:,k)
  enddo
  ! phase velocity
  c1 = 0d0; c1W = 0d0
  do k=1,nz-1 
    c1(:,:)=c1(:,:)+sqrt(Nsqr_lim(:,:,k))*dzw(k)*maskW(:,:,k)/pi
  enddo
  do k=1,nz  
    c1W(:,:)=c1W(:,:)+sqrt(Nsqrt(:,:,k))*dzt(k)*maskT(:,:,k)/pi
  enddo  
  ! Rossby radius 
  R1 = maskT(:,:,nz)*c1/max(rossmix2_fmin,abs(coriolis_t) )  
  ! calculate xi=int_(-h)^z N dz, xi is on T grid
  xi=0; xiW=0
  do k=2,nz   
   xi(:,:,k)  = xi(:,:,k-1)*maskT(:,:,k-1)  + sqrt(Nsqr_lim(:,:,k-1))*dzw(k-1)*maskW(:,:,k-1) 
   xiW(:,:,k) = xiW(:,:,k-1)*maskW(:,:,k-1) + sqrt(Nsqrt(:,:,k))*dzt(k)*maskT(:,:,k) 
  enddo
  ! calculate phi_n    \sim cos(xi/c_n)*N^0.5 and  (f/N)dphi_n/dz  \sim sin(xi/c_n)*N^0.5
  do k=1,nz   
    phi1(:,:,k)  = cos(   xi(:,:,k)/max(1d-22,c1) )*Nsqrt(:,:,k)**0.25*maskT(:,:,k)
    phi1W(:,:,k) = cos(  xiW(:,:,k)/max(1d-22,c1W) )*Nsqr_lim(:,:,k)**0.25*maskW(:,:,k)
    phiz1(:,:,k)  = sin(   xi(:,:,k)/max(1d-22,c1(:,:)) )*Nsqrt(:,:,k)**0.25*maskT(:,:,k)
    phiz1W(:,:,k) = sin(  xiW(:,:,k)/max(1d-22,c1W(:,:)) )*Nsqr_lim(:,:,k)**0.25 *maskW(:,:,k)
  enddo
    
  ! normalisation of phi_n  and  (f/N)dphi_n/dz 
  ! int_(-h)^0 dz phi_n^2  = 1 
  norm=0
  do k=1,nz   
    norm = norm + phi1(:,:,k)**2*dzt(k)*maskT(:,:,k)
  enddo
  do k=1,nz   
   where( norm/=0) phi1(:,:,k) = phi1(:,:,k)/norm**0.5
  enddo
  ! int_(-h)^0  ( (f/N)dphi_n/dz)^2 dz  = 1/R^2_n
  norm=0
  do k=1,nz   
    norm = norm + phiz1(:,:,k)**2*dzt(k)*maskT(:,:,k)
  enddo
  do k=1,nz   
     where( norm/=0) phiz1(:,:,k) = phiz1(:,:,k)/norm**0.5/max(1d-12,R1(:,:))
  enddo
 
  ! also on W grid  
  norm=0
  do k=1,nz   
   norm = norm+ phi1W(:,:,k)**2*dzw(k)*maskW(:,:,k)
  enddo
  do k=1,nz   
   where( norm/=0) phi1W(:,:,k) = phi1W(:,:,k)/norm**0.5
  enddo
  norm=0
  do k=1,nz   
   norm = norm + phiz1W(:,:,k)**2*dzw(k)*maskW(:,:,k)
  enddo
  do k=1,nz   
   where( norm/=0) phiz1W(:,:,k) = phiz1W(:,:,k)/norm**0.5/max(1d-12,R1(:,:))
  enddo    
  
end subroutine rossmix2_vertical_mode 





subroutine rossmix2_geo_velocity
 ! calculate geostrophic velocities
 use main_module   
 use rossmix2_module   
 implicit none
 integer :: i,j,k,itts_loc
 real*8 :: pg(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz),fxa,f0,fxb
 real*8 :: fpx(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx)
 real*8 :: fpy(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx)
 real*8 :: forc(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx)
 real*8,save,allocatable :: spress(:,:,:)
 logical, save :: first = .true.
 
 if (first) then
  !hydrostatic pressure
  fxa = grav/rho_0
  p_hydro(:,:,nz) = 0.5*rho(:,:,nz,tau)*fxa*dzw(nz)*maskT(:,:,nz)
  do k=nz-1,1,-1
   p_hydro(:,:,k)= maskT(:,:,k)*(p_hydro(:,:,k+1)+ 0.5*(rho(:,:,k+1,tau)+rho(:,:,k,tau))*fxa*dzw(k))
  enddo
  allocate( spress(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,3) ); spress=0
 endif
 
 
 if (enable_streamfunction) then  
  if (my_pe==0) print*,' diagnosing surface pressure for geostrophic velocities'
  ! forcing for surface pressure 
  fpx=0.;fpy=0.
  do k=1,nz
    do j=js_pe,je_pe
      do i=is_pe,ie_pe
       fxa = u(i,j,k,tau)+dt_mom*( du_mix(i,j,k)+ (1.5+AB_eps)*du(i,j,k,tau)  &
                                 - (0.5+AB_eps)*du(i,j,k,taum1) )*maskU(i,j,k)       
       fxb = v(i,j,k,tau)+dt_mom*( dv_mix(i,j,k)+ (1.5+AB_eps)*dv(i,j,k,tau)  &
                                 - (0.5+AB_eps)*dv(i,j,k,taum1) )*maskV(i,j,k)                                 
       fpx(i,j)=fpx(i,j) + fxa*maskU(i,j,k)*dzt(k)/dt_mom                          
       fpy(i,j)=fpy(i,j) + fxb*maskV(i,j,k)*dzt(k)/dt_mom
      enddo
     enddo
  enddo
  call border_exchg_xy(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,fpx) 
  call setcyclic_xy   (is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,fpx)
  call border_exchg_xy(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,fpy)
  call setcyclic_xy   (is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,fpy)
  ! forc = 1/cos (u_x + (cos v)_y )A
  do j=js_pe,je_pe
   do i=is_pe,ie_pe
    forc(i,j)=(fpx(i,j)-fpx(i-1,j))/(cost(j)*dxt(i))+(cosu(j)*fpy(i,j)-cosu(j-1)*fpy(i,j-1))/(cost(j)*dyt(j))
   enddo
  enddo
  if (enable_free_surface) then
   do j=js_pe,je_pe
    do i=is_pe,ie_pe
     forc(i,j) = forc(i,j) - spress(i,j,tau)/(grav*dt_mom**2)*maskT(i,j,nz)
    enddo
   enddo
  endif 
  !solve for surface pressure
  fpx = psi(:,:,taup1); fxa = congr_epsilon 
  psi(:,:,taup1) = 2*spress(:,:,tau) - spress(:,:,taum1) ! first guess
  congr_epsilon = rossmix2_congr_epsilon
  call congrad_surf_press(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,forc,itts_loc)
  call border_exchg_xy(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,psi(:,:,taup1)); 
  call setcyclic_xy(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,psi(:,:,taup1))
  if (my_pe==0) print*,' itts = ',itts_loc
  spress(:,:,taup1) = psi(:,:,taup1)
  psi(:,:,taup1) = fpx; congr_epsilon = fxa
 else
  spress(:,:,taup1) = psi(:,:,taup1)
 endif
 
  
 ! geostrophic velocities from pressure gradients
 do k=1,nz
  pg(:,:,k) = p_hydro(:,:,k) + spress(:,:,taup1)
 enddo
 do j=js_pe,je_pe 
   do i=is_pe,ie_pe
    f0 = sign(max(rossmix2_fmin,abs(coriolis_t(i,j))),coriolis_t(i,j))
    ug(i,j,:) = 0.5*(-( pg(i,j+1,:)-pg(i,j  ,:))/dyu(j  )*maskV(i,j  ,:) &
                     -( pg(i,j  ,:)-pg(i,j-1,:))/dyu(j-1)*maskV(i,j-1,:) ) /f0*maskT(i,j,:)
    vg(i,j,:) = 0.5*(( pg(i+1,j,:)-pg(i  ,j,:))/(dxu(i  )*cost(j))*maskU(i  ,j,:) & 
                   + ( pg(i  ,j,:)-pg(i-1,j,:))/(dxu(i-1)*cost(j))*maskU(i-1,j,:) )/f0*maskT(i,j,:)
   enddo
 enddo 
 first = .false.
end subroutine rossmix2_geo_velocity





subroutine rossmix2_appr_wavenumber
 use main_module   
 use rossmix2_module   
 implicit none
 integer :: i,j,k,ks
 real*8 :: f0,h0,um,vm,uvm,uum,vvm,theta,uk(nz)
 complex*16 :: c0,ps1(nz),pp(nz),q
 
 do j=js_pe,je_pe
  do i=is_pe,ie_pe
   ks = kbot(i,j)
   if (ks>0) then
     f0 = sign( max(rossmix2_fmin, abs(coriolis_t(i,j)) ), coriolis_t(i,j))
     h0 = max(1d-12,ht(i,j))

     ! direction of approximate growth rate maximum    
     um = sum(ug(i,j,:)*dzt)/h0
     vm = sum(vg(i,j,:)*dzt)/h0
     uvm= sum(ug(i,j,:)*vg(i,j,:)*dzt)/h0
     uum = sum(ug(i,j,:)**2*dzt)/h0
     vvm = sum(vg(i,j,:)**2*dzt)/h0
     theta = 0.5*atan2(2*(uvm-um*vm), uum-um**2-vvm+vm**2)
     uk=ug(i,j,:)*cos(theta)+vg(i,j,:)*sin(theta)
 
     ! magnitude of wave number of approximate minimum
     c0 = um**2 - uum
     c0 = um + sqrt(c0)
     ps1=0. 
     do k=ks+1,nz
      ps1(k) = ps1(k-1)+dzt(k)*(c0-uk(k))**2 
     enddo
     pp=0
     do k=ks+1,nz
      pp(k) = pp(k-1)+dzt(k)*( Nsqr_lim(i,j,k)/(f0**2*(c0-uk(k))**2)*ps1(k) )
     enddo
     q = 2*c0*sum(Uk*pp*dzt)/h0 -c0**2*sum(pp*dzt)/h0 - sum(Uk**2*pp*dzt)/h0
     q = q/( 2*(c0-sum(Uk*dzt)/h0))
    
     kh_appr(i,j)=  sqrt(abs( imag(c0)/(3*imag(q))  ))
     kh_appr(i,j) = min(rossmix2_kh_appr_max/max(1d-12,R1(i,j)), kh_appr(i,j) )  ! limit wavenumber
     kh_appr(i,j) = max(rossmix2_kh_appr_min/max(1d-12,R1(i,j)), kh_appr(i,j) )  ! limit wavenumber
   endif 
  enddo
 enddo 
end subroutine rossmix2_appr_wavenumber





subroutine rossmix2_eigenvalue_simple
 ! Rossmix1 vertical structure functions
 use main_module   
 use rossmix2_module   
 implicit none
 integer :: i,j,k,ks,n
 real*8 :: f0,fxa
    
 if (my_pe==0) print*,'solving simple local eigenvalue problems at itt=',itt 
 ub_lim = 0d0
 uu_lim = 0d0
 do j=js_pe,je_pe
  do i=is_pe,ie_pe
   ks = kbot(i,j)
   if (ks>0 .and. ht(i,j)>rossmix2_mindepth) then
     f0 = sign( max(rossmix2_fmin, abs(coriolis_t(i,j)) ), coriolis_t(i,j))    
     do k=ks,nz-1   
      fxa = phi1W(i,j,k)**2*(rossmix2_A8*R1(i,j)) + phiz1W(i,j,k)**2*(rossmix2_A7*R1(i,j)**3)
      ub(i,j,k,:)= 2*fxa*sqrt(Nsqr_lim(i,j,k))/R1(i,j)/f0*maskW(i,j,k)
     enddo        
     do n=2,nphi-1
      ub(i,j,ks:nz-1,n) = ub(i,j,ks:nz-1,n) - ub(i,j,ks,n)
     enddo 
     ub(i,j,nz,:)= 0d0 
     do n=2,nphi-1
      uu(i,j,ks:nz,n)= 2*rossmix2_A2*phi1(i,j,ks:nz)**2
     enddo 
 
      ! diagnose limiter for structure functions 
     do n=2,nphi-1 
      do k=1,nz
       if ( abs(ub(i,j,k,n)) > rossmix2_ub_limiter*Nsqr_lim(i,j,k) ) &
                 ub_lim(i,j,k) = ub_lim(i,j,k) + dphit(n)
       if ( abs(uu(i,j,k,n)) > rossmix2_uu_limiter )  uu_lim(i,j,k) = uu_lim(i,j,k) + dphit(n)          
      enddo
     enddo 
     
     ! limit structure functions     
     do k=1,nz
       ub(i,j,k,:)=sign( min(rossmix2_ub_limiter*Nsqr_lim(i,j,k), abs(ub(i,j,k,:)) ), ub(i,j,k,:)) 
       uu(i,j,k,:)=sign( min(rossmix2_uu_limiter, abs(uu(i,j,k,:)) ), uu(i,j,k,:)) 
     enddo       
   endif   
  enddo
 enddo
 
 do k=1,nz
  call border_exchg_xyp(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nphi,ub(:,:,k,:)) 
  call setcyclic_xyp(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nphi,ub(:,:,k,:))
  call border_exchg_xyp(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nphi,uu(:,:,k,:)) 
  call setcyclic_xyp(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nphi,uu(:,:,k,:))
 enddo

 if (enable_rossmix2_smooth_mom_fluxes) call rossmix2_smooth(uu,2) 
 
 if (my_pe==0) print*,'done solving local eigenvalue problems'
end subroutine rossmix2_eigenvalue_simple







subroutine rossmix2_eigenvalue_approx
 ! calculate vertical eigenfunction and eigenvalue approximately
 use main_module   
 use rossmix2_module   
 implicit none
 integer :: i,j,k,ks,n
 real*8 :: kx,ky,norm,kh,kh2,uk(nz),betas,f0,h0,bar_u,fxa
 complex*16 :: c0,c1,h1s(nz),h1(nz),h2(nz),cmax, phi(nz),phiz(nz),h2s(nz),c2  
  
 call rossmix2_geo_velocity
  
 if (enable_rossmix2_appr_wavenumber) then
  call rossmix2_appr_wavenumber
 else
  kh_appr(:,:) = (1.4/pi)/max(1d-12,R1)
 endif
 
 if (my_pe==0) print*,'solving approx. local eigenvalue problems at itt=',itt 
 ub_lim = 0d0
 uu_lim = 0d0
 do j=js_pe,je_pe
  do i=is_pe,ie_pe
   ks = kbot(i,j)
   if (ks>0 .and. ht(i,j)>rossmix2_mindepth) then
    do n=2,nphi-1
     f0 = sign( max(rossmix2_fmin, abs(coriolis_t(i,j)) ), coriolis_t(i,j))
     h0 = max(1d-12,ht(i,j))

     kx = kh_appr(i,j)*cos(phit(n))
     ky = kh_appr(i,j)*sin(phit(n))
   
     kh2 = kx**2+ky**2
     kh = sqrt(kh2)
     uk = (ug(i,j,:)*kx+vg(i,j,:)*ky)/kh
     betas = beta(i,j) * kx/kh

     ! approx. eigenvalue
     bar_u = sum(uk*dzt)/h0
     c0 = bar_u**2 - sum(uk**2*dzt)/h0 + (betas/(2*kh2))**2 
     c0 = bar_u - betas/(2*kh2) + sqrt(c0) 
     h1s=0;h1=0    
     do k=ks+1,nz
      h1s(k)=h1s(k-1) + dzt(k)*( (c0-Uk(k))**2*kh2 + betas* (c0-Uk(k)))  
     enddo
     do k=ks+1,nz
      h1(k)=h1(k-1) + dzt(k)*( h1s(k)*Nsqr_lim(i,j,k)/( f0**2*(c0-Uk(k))**2 ) )  
     enddo
     h1 = h1*maskT(i,j,:)
     c1 = c0*(c0+betas/kh2)*sum(h1*dzt)/h0 &
            -(2*c0+betas/kh2)*sum(h1*Uk*dzt)/h0 + sum(h1*Uk**2*dzt)/h0
     c1 =-c1/(betas/kh2+2*(c0-bar_u) ) 
     cmax = c0 + c1
     
if (.false.) then     
     h2s=0;h2=0    
     do k=ks+1,nz
      h2s(k)=h2s(k-1) + dzt(k)*H1(k)*( (c0+c1-Uk(k))**2*kh2 + betas* (c0+c1-Uk(k)))  
     enddo
     do k=ks+1,nz
      h2(k)=h2(k-1) + dzt(k)*( h2s(k)*Nsqr_lim(i,j,k)/( f0**2*(c0+c1-Uk(k))**2 ) )  
     enddo
     h2 = h2*maskT(i,j,:)
     
     c2 =   (sum(h2*Uk*dzt)/h0 - sum(h1*dzt)/h0*c1)*(2*c0+betas/kh2 ) &
          - sum(h2*dzt)/h0*c0*(betas/kh2+c0) -c1**2 +2*c1*sum(h1*Uk*dzt)/h0 - sum(h2*Uk**2*dzt)/h0 
     c2 = c2/(betas/kh2+2*(c0-bar_u) )
     
     cmax = c0 + c1 + c2
endif     
     
     om_max(i,j,n) = cmax*kh 
 
     if (abs(imag(om_max(i,j,n)))>0) then
      ! approximate eigenvector
      H1(:)=0.
      do k=ks+1,nz
        H1(k)=H1(k-1)+dzt(k)*(uk(k)-cmax)*(kh2*(uk(k)-cmax)-betas) 
      enddo
      do k=nz+1,nz
       H1(k) = H1(k-1)+dzt(k)*Nsqr_lim(i,j,k)/f0**2/(uk(k)-cmax)**2*H1(k)
      enddo
      H2(:)=0.0
      do k=ks+1,nz
        H2(k)=H2(k-1)+dzt(k)*H1(k)*(cmax-uk(k))*(kh2*(uk(k)-cmax)-betas) 
      enddo
      do k=ks+1,nz
       H2(k) = H2(k-1)+dzt(k)*Nsqr_lim(i,j,k)/f0**2/(uk(k)-cmax)**2*H2(k)
      enddo      
      phi = (1+H1+H2)*(uk-cmax)
      
      ! vertical derivative and normalisation
      do k=ks,nz-1
       phiz(k) = (phi(k+1)-phi(k))/dzw(k)*maskW(i,j,k) 
      enddo
      phiz(nz) =  0.0  
      norm=0
      do k=ks,nz   
       ! here coriolis_t should be f0
       norm = norm + (kx**2+ky**2)*abs(phi(k))**2*dzt(k)*maskT(i,j,k) &
                   + coriolis_t(i,j)**2/Nsqr_lim(i,j,k)*abs(phiz(k))**2*dzw(k)*maskW(i,j,k)
      enddo     
      if (norm/=0) phi = phi/norm**0.5*(kx**2+ky**2+1/R1(i,j)**2)**0.5
      if (norm/=0) phiz = phiz/norm**0.5*(kx**2+ky**2+1/R1(i,j)**2)**0.5
      
      ! structure function 
      do k=ks,nz-1   
       ub(i,j,k,n)= 2*rossmix2_A1*R1(i,j)* &
                    real(cmplx(0,1)*0.5*(phi(k)+phi(k+1))*conjg(phiz(k)))*maskW(i,j,k)
      enddo     
      ub(i,j,ks:nz-1,n) = ub(i,j,ks:nz-1,n) - ub(i,j,ks,n)
      ub(i,j,nz,n)= 0d0     
      uu(i,j,ks:nz,n)= 2*rossmix2_A2*abs(phi(ks:nz))**2  
     else      
      ! fall back to Rossmix1 vertical structure function
      do k=ks,nz-1         
       fxa = phi1W(i,j,k)**2*(rossmix2_A8*R1(i,j)) + phiz1W(i,j,k)**2*(rossmix2_A7*R1(i,j)**3)       
       ub(i,j,k,n)= 2*fxa*sqrt(Nsqr_lim(i,j,k))/R1(i,j)/f0*maskW(i,j,k)      
      enddo        
      ub(i,j,ks:nz-1,n) = ub(i,j,ks:nz-1,n) - ub(i,j,ks,n)
      ub(i,j,nz,n)= 0d0 
      uu(i,j,ks:nz,n)= 2*rossmix2_A2*phi1(i,j,ks:nz)**2
     endif 
 
     ! diagnose limiter 
     do k=1,nz
       if ( abs(ub(i,j,k,n)) > rossmix2_ub_limiter*Nsqr_lim(i,j,k) ) &
                 ub_lim(i,j,k) = ub_lim(i,j,k) + dphit(n)
       if ( abs(uu(i,j,k,n)) > rossmix2_uu_limiter )  uu_lim(i,j,k) = uu_lim(i,j,k) + dphit(n)          
     enddo
     
     ! B = int dphi E ub/N^2 f n , ub = 2 A_ 1 real(i*phi*conjg(phiz))
     ! B =  int dphi 2 E f A_ 1 real(i*phi*conjg(phiz))  n /N^2
     ! limit structure functions
     do k=1,nz
       ub(i,j,k,n)=sign( min(rossmix2_ub_limiter*Nsqr_lim(i,j,k), abs(ub(i,j,k,n)) ), ub(i,j,k,n)) 
       uu(i,j,k,n)=sign( min(rossmix2_uu_limiter, abs(uu(i,j,k,n)) ), uu(i,j,k,n)) 
     enddo       
          
    enddo
   endif   
  enddo
 enddo
 
 do k=1,nz
  call border_exchg_xyp(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nphi,ub(:,:,k,:)) 
  call setcyclic_xyp(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nphi,ub(:,:,k,:))
  call border_exchg_xyp(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nphi,uu(:,:,k,:)) 
  call setcyclic_xyp(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nphi,uu(:,:,k,:))
 enddo
 if (enable_rossmix2_smooth_mom_fluxes) call rossmix2_smooth(uu,2) 

 if (my_pe==0) print*,'done solving local eigenvalue problems'
end subroutine rossmix2_eigenvalue_approx





subroutine rossmix2_eigenvalue_lapack
 ! solve quasi geostrophic vertical eigen value problem after
 ! Smith (2007) The Geography of Linear Baroclinic 
 ! Instability in Earth's Oceans, J. Mar. Res. 
 use main_module   
 use rossmix2_module   
 implicit none
 complex*16 :: phi(nz),phiz(nz)
 integer :: i,j,k,ks,n,k_max,info,nn
 real*8 :: G_(nz,nz),Qx(nz),Qy(nz),B_(nz,nz),A_(nz,nz),kx,ky,norm,f0,fxa
 real*8 :: alphar(nz),alphai(nz),beta_loc(nz),vl(nz,nz),vr(nz,nz)
 real*8, allocatable, save :: work(:)
 integer, save :: lwork
 logical, save :: first = .true.
 
 ! determine optimal work space for lapack solver
 if (first) then
  allocate(work(10))
  call dggev('N','V',nz,A_,nz,B_,nz,alphar,alphai,beta_loc,vl,nz,vr,nz,work,-1,info)
  lwork = int(work(1)) + 100
  call global_max_int(lwork)
  deallocate(work)
  allocate(work(lwork)); work=0.0  
 endif

 call rossmix2_geo_velocity
 
 if (enable_rossmix2_appr_wavenumber) then
  call rossmix2_appr_wavenumber
 else
  kh_appr(:,:) = (1.4/pi)/max(1d-12,R1)
 endif
 
 if (my_pe==0) print*,'solving local eigenvalue problems at itt=',itt 
 ub_lim = 0d0
 uu_lim = 0d0
 do j=js_pe,je_pe
  do i=is_pe,ie_pe
   ks = kbot(i,j)
   if (ks>0 .and. ht(i,j)>rossmix2_mindepth) then
    f0 = sign( max(rossmix2_fmin, abs(coriolis_t(i,j)) ), coriolis_t(i,j))
    do n=2,nphi-1

     kx = kh_appr(i,j)*cos(phit(n))
     ky = kh_appr(i,j)*sin(phit(n))
          
     ! vertical differencing operator     
     G_ = 0.         
     G_(ks,ks)   = -f0**2*( 1.0/(Nsqr_lim(i,j,ks)*dzw(ks)))/dzw(ks)*maskW(i,j,ks)
     G_(ks,ks+1) = -f0**2*(-1.0/(Nsqr_lim(i,j,ks)*dzw(ks)))/dzw(ks)*maskW(i,j,ks)
     do k=ks+1,nz-1
      G_(k,k-1) = -f0**2*(-1./(Nsqr_lim(i,j,k-1)*dzw(k-1)))/dzw(k)  *maskW(i,j,k)  
      G_(k,k)   = -f0**2*( 1./(Nsqr_lim(i,j,k-1)*dzw(k-1))+1./(Nsqr_lim(i,j,k)*dzw(k)))/dzw(k)   *maskW(i,j,k) 
      G_(k,k+1) = -f0**2*(-1./(Nsqr_lim(i,j,k)*dzw(k)))/dzw(k) *maskW(i,j,k)
     enddo  
     G_(nz,nz)   = -f0**2*( 1./(Nsqr_lim(i,j,k-1)*dzw(k-1)))/dzw(nz) *maskW(i,j,nz)
     G_(nz,nz-1) = -f0**2*(-1./(Nsqr_lim(i,j,k-1)*dzw(k-1)))/dzw(nz) *maskW(i,j,nz)
     
        
     ! background PV gradient   
     Qy = beta(i,j)*maskT(i,j,:)-matmul(G_,ug(i,j,:)) 
     Qx = matmul(G_,vg(i,j,:))
     ! topography
     !Qx[-1]+=f0*hx/dz
     !Qy[-1]+=f0*hy/dz
     
     B_ = G_  
     do k=1,nz
      B_(k,k) = B_(k,k) - (kx**2+ky**2)*maskT(i,j,k)
     enddo
   
     A_ = 0.
     do k=1,nz
      A_(k,k) = (kx*ug(i,j,k) + ky*vg(i,j,k))*maskT(i,j,k)
     enddo
     A_ = matmul(A_,B_)
     do k=1,nz
      A_(k,k) = A_(k,k) + (kx*Qy(k) - ky*Qx(k))*maskT(i,j,k)
     enddo
     ! omega B_ij v_j = A_ij v_j     
     nn=nz-ks+1
     call dggev('N','V',nn,A_(ks:nz,ks:nz),nn,B_(ks:nz,ks:nz),nn,alphar,alphai,beta_loc, &
                           vl,nz,vr,nz,work,lwork,info)
                      
     !if (info/=0) print*,'WARNING: info=',info,' in rossmix2_eigenvalue at i,j=',i,j
     alphar = alphar/beta_loc
     alphai = alphai/beta_loc
                                         
     ! pick maximum growth rate     
     k_max= maxloc(alphai(1:nn),1)
     om_max(i,j,n) = cmplx(alphar(k_max),alphai(k_max),kind=8)  
          
     if (k_max<nn.and.abs(alphai(k_max))>0.and.info==0) then  
       do k=1,nn
        phi(k+ks-1) = cmplx(VR(k,k_max),VR(k,k_max+1),kind=8)        
       enddo 
       phi = phi*maskT(i,j,:)
       ! vertical derivative and normalisation
       do k=ks,nz-1
        phiz(k) = (phi(k+1)-phi(k))/dzw(k)*maskW(i,j,k) 
       enddo
       phiz(nz) =  0.0    
       norm=0
       do k=ks,nz   
        norm = norm +  (kx**2+ky**2)*abs(phi(k))**2*dzt(k)*maskT(i,j,k) 
       enddo
       do k=ks,nz-1
        norm = norm + f0**2/Nsqr_lim(i,j,k)*abs(phiz(k))**2*dzw(k)*maskW(i,j,k)
       enddo     
       do k=ks,nz   
        if (norm/=0) phi(k) = phi(k)/norm**0.5*(kx**2+ky**2+1/R1(i,j)**2)**0.5
        if (norm/=0) phiz(k) = phiz(k)/norm**0.5*(kx**2+ky**2+1/R1(i,j)**2)**0.5
       enddo

       ! structure function 
       do k=ks,nz-1   
         ub(i,j,k,n)= 2*rossmix2_A1*R1(i,j)* &
                      real(cmplx(0,1)*0.5*(phi(k)+phi(k+1))*conjg(phiz(k)))*maskW(i,j,k)            
       enddo     
       ub(i,j,ks:nz-1,n) = ub(i,j,ks:nz-1,n) - ub(i,j,ks,n)
       ub(i,j,nz,n)= 0d0
       uu(i,j,ks:nz,n)= 2*rossmix2_A2*abs(phi(ks:nz))**2 
     
     else
      ! fall back to Rossmix1 vertical structure function
      do k=ks,nz-1         
       fxa = phi1W(i,j,k)**2*(rossmix2_A8*R1(i,j)) + phiz1W(i,j,k)**2*(rossmix2_A7*R1(i,j)**3)       
       ub(i,j,k,n)= 2*fxa*sqrt(Nsqr_lim(i,j,k))/R1(i,j)/f0*maskW(i,j,k)      
      enddo        
      ub(i,j,ks:nz-1,n) = ub(i,j,ks:nz-1,n) - ub(i,j,ks,n)
      ub(i,j,nz,n)= 0d0 
      uu(i,j,ks:nz,n)= 2*rossmix2_A2*phi1(i,j,ks:nz)**2       
     endif

     ! diagnose limiter for structure functions 
     do k=1,nz
       if ( abs(ub(i,j,k,n)) > rossmix2_ub_limiter*Nsqr_lim(i,j,k) ) &
                 ub_lim(i,j,k) = ub_lim(i,j,k) + dphit(n)
       if ( abs(uu(i,j,k,n)) > rossmix2_uu_limiter )  uu_lim(i,j,k) = uu_lim(i,j,k) + dphit(n)          
     enddo
 
     ! limit structure functions
     do k=1,nz
      ub(i,j,k,n)=sign( min(rossmix2_ub_limiter*Nsqr_lim(i,j,k), abs(ub(i,j,k,n)) ), ub(i,j,k,n)) 
      uu(i,j,k,n)=sign( min(rossmix2_uu_limiter, abs(uu(i,j,k,n)) ), uu(i,j,k,n)) 
     enddo       

    enddo
   endif   
  enddo
 enddo
 
 do k=1,nz
  call border_exchg_xyp(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nphi,ub(:,:,k,:)) 
  call setcyclic_xyp(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nphi,ub(:,:,k,:))
  call border_exchg_xyp(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nphi,uu(:,:,k,:)) 
  call setcyclic_xyp(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nphi,uu(:,:,k,:))
 enddo
 if (enable_rossmix2_smooth_mom_fluxes) call rossmix2_smooth(uu,2)

 if (my_pe==0) print*,'done solving local eigenvalue problems'
 first = .false.
end subroutine rossmix2_eigenvalue_lapack




subroutine rossmix2_smooth(u_,ij)
 ! smooth u_ over ijxij grid points
 use main_module   
 use rossmix2_module   
 implicit none
 integer :: n,ii,jj,i,j,k,ij
 real*8 :: u_(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz,nphi)
 real*8 :: bloc(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz)
 real*8 :: aloc(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz)

 do n=2,nphi-1
  bloc=0d0; aloc=0d0
  do jj=-ij,ij
   do ii=-ij,ij
    do j=js_pe,je_pe
     do i=is_pe,ie_pe
      bloc(i,j,:) = bloc(i,j,:) + u_(i+ii,j+jj,:,n)*maskT(i+ii,j+jj,:)
      aloc(i,j,:) = aloc(i,j,:) + maskT(i+ii,j+jj,:)
     enddo
    enddo    
   enddo
  enddo
  u_(:,:,:,n) = bloc/max(1d-12,aloc)
 enddo 
 do k=1,nz
  call border_exchg_xyp(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nphi,u_(:,:,k,:)) 
  call setcyclic_xyp(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nphi,u_(:,:,k,:))
 enddo
end subroutine rossmix2_smooth


