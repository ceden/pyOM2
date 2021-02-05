



subroutine rossmix_main
 use main_module   
 use rossmix_module   
 implicit none
 integer :: m
 logical, save :: init = .false.

 if (.not. enable_rossmix) return

 if ( .not. init .or. mod(itt,int(rossmix_calc_cg_int/dt_tracer))  == 0) then

  init = .true.
  !---------------------------------------------------------------------------------
  ! calc group velocity etc only rarely
  !---------------------------------------------------------------------------------

  call calc_rossmix_vertical_struct_fct
  call set_rossmix_parameter

  do m=1,nmodes
   call rossmix_group_velocity(m)
   call rossmix_check_cfl(m)
   if (enable_rossmix_mean_flow_interaction) call rossmix_growth_rate(m)
  enddo

 endif

 if (enable_rossmix_mean_flow_interaction) then
    call rossmix_eddy_streamfunction
    call rossmix_eddy_velocity
    if (enable_rossmix_lateral_stress)  call rossmix_lateral_stress
 endif

 call integrate_rossmix

end subroutine rossmix_main






subroutine rossmix_initialize
  use main_module   
  use rossmix_module   
  implicit none
  integer :: i,j,k

  if (.not. enable_rossmix) return

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

  call reflect_rossmix_initialize

end subroutine rossmix_initialize






subroutine set_rossmix_parameter
 use main_module   
 use rossmix_module   
 implicit none
 integer :: i,j,m
 real*8 :: rRn,rRx,fxc,fxc2,afxc,lfxc2
 real*8 :: Rx(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx)
 real*8 :: A0_l,A0_s,A1_l,A1_s, A2_l,A2_s
 real*8 :: A3_l,A3_s,A4_l,A4_s, A6_l,A6_s
 real*8 :: A7_l,A7_s,A8_l,A8_s


 do m=1,nmodes
  Rx(:,:) = 2*Rn(:,:,2) 

  do j=js_pe-onx,je_pe+onx
   do i=is_pe-onx,ie_pe+onx
 
    rRn = 1./max(1d-12,Rn(i,j,m))
    rRx = 1./max(1d-12,Rx(i,j))
    fxc  = max(1d-18,Rx(i,j)*rRn )
    fxc2 = fxc**2 
    afxc = atan(fxc)
    lfxc2 = log(fxc2)

    km_l(i,j,m) = log(fxc2+1)/(2*afxc)*rRx

    A0_s = (pi/2-afxc)*rRx  
    A0_l = afxc*rRx   
    A1_l = log(4./(fxc2+1)**2)/( 4*(1-fxc2) )
    A1_s = ( 2*log(fxc2)+log(4./(fxc2+1)**2) )/( 4*(fxc2-1) )
    A2_s = rRx*((pi/4)*fxc+afxc -pi/2)/ (fxc2-1) 
    A2_l = rRx*((pi/4)*fxc-afxc )     / (fxc2-1)
    A3_s = Rn(i,j,m)*( (pi-2)*fxc2 +8*fxc*afxc +pi+2 -4*pi*fxc )/( 8*(fxc2-1)**2 )
    A3_l = Rn(i,j,m)*( (pi+2)*fxc2 -8*fxc*afxc +pi-2           )/( 8*(fxc2-1)**2 )
    A4_s = Rn(i,j,m)*(pi+4*fxc*afxc -2*pi*fxc)/( 4*(1-fxc2) )
    A4_l = Rn(i,j,m)*(pi-4*fxc*afxc          )/( 4*(1-fxc2) )
    A6_s = Rn(i,j,m)**3*(fxc2*(4*pi*fxc -8*fxc*afxc-3*pi+2) +pi-2)/( 8*(fxc2-1)**2 )
    A6_l = Rn(i,j,m)**3*(fxc2*(          8*fxc*afxc-3*pi-2) +pi+2)/( 8*(fxc2-1)**2 )

    A7_l = Rn(i,j,m)**2*(1-fxc2*(            -2*log(fxc2+1)+1+log(4.) ) )/( 4*(fxc2-1)**2 )
    A7_s = Rn(i,j,m)**2*(1+fxc2*( 2*lfxc2-2*log(fxc2+1)-1+log(4.) ) )/( 4*(fxc2-1)**2 )
    A8_l = (fxc2            -2*log(fxc2+1) -1 +log(4.) )/( 4*(fxc2-1)**2 )
    A8_s = (fxc2-2*lfxc2+2*log(fxc2+1) -1 -log(4.) )/( 4*(fxc2-1)**2 )

    cfac1_l(i,j,m) =  A4_l/A0_l
    cfac1_s(i,j,m) =  A4_s/A0_s
    cfac2_l(i,j,m) =  A3_l/A0_l
    cfac2_s(i,j,m) =  A3_s/A0_s
   
    if (m==1) then ! switch off barotropic long waves
     phi_ga_s(i,j,:,m) = phinW(i,j,:,m)**2*A1_s/A0_s
    else
     phi_ga_l(i,j,:,m) = (phinW(i,j,:,m)**2*A8_l + phiznW(i,j,:,m)**2*A7_l )/A0_l
     phi_ga_s(i,j,:,m) = (phinW(i,j,:,m)**2*A8_s + phiznW(i,j,:,m)**2*A7_s )/A0_s
     phi_nu_l(i,j,m)   = A2_l/A0_l
    endif
    phi_nu_s(i,j,m)   = A2_s/A0_s

    ray_fric_l(i,j,m) = A2_l/A0_l*rossmix_rayleigh_friction
    ray_fric_s(i,j,m) = A2_s/A0_s*rossmix_rayleigh_friction

   enddo
  enddo
 enddo

end subroutine set_rossmix_parameter








subroutine calc_rossmix_vertical_struct_fct
  ! calculate vertical structure function
  ! Nsqr and frequencies must be set before this routine
  use main_module   
  use rossmix_module   
  implicit none
  integer :: k,n
  real*8 :: norm(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx)
  real*8 :: xi(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz)
  real*8 :: xiW(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz)
  real*8 :: Nsqr_lim(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz)
  real*8 :: Nsqrt(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz)
  real*8 :: cn(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nmodes)
  real*8 :: cnW(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nmodes)

  ! limit and interpolate Nsqr
  Nsqr_lim = max(rossmix_N2min,Nsqr(:,:,:,tau))
  Nsqrt(:,:,1) = Nsqr_lim(:,:,1)
  do k=2,nz
    Nsqrt(:,:,k) = (Nsqr_lim(:,:,k)*maskW(:,:,k)+Nsqr_lim(:,:,k-1)*maskW(:,:,k-1)) / &
                                (maskW(:,:,k)+maskW(:,:,k-1)+1d-22)*maskT(:,:,k)
  enddo


  ! phase velocity
  cn = 0.0 ; cnW = 0.0
  cn(:,:,1)  = sqrt(9.81*max(rossmix_hmin,ht)) !barotropic wave speed 
  cnW(:,:,1) = cn(:,:,1) 
  do n=2,nmodes
   do k=1,nz-1 
    cn(:,:,n)=cn(:,:,n)+sqrt(Nsqr_lim(:,:,k))*dzw(k)*maskW(:,:,k)/( (n-1)*pi)
   enddo
   do k=1,nz  
    cnW(:,:,n)=cnW(:,:,n)+sqrt(Nsqrt(:,:,k))*dzt(k)*maskT(:,:,k)/( (n-1)*pi)
   enddo
  enddo

  ! Rossby radius
  do n=1,nmodes
   Rn(:,:,n) = maskT(:,:,nz)*cn(:,:,n)/max(rossmix_fmin,abs( coriolis_t(:,:) ) )
  enddo


  ! calculate xi=int_(-h)^z N dz, xi is on T grid
  xi=0; xiW=0
  do k=2,nz   
   xi(:,:,k)  = xi(:,:,k-1)*maskT(:,:,k-1)  + sqrt(Nsqr_lim(:,:,k-1))*dzw(k-1)*maskW(:,:,k-1) 
   xiW(:,:,k) = xiW(:,:,k-1)*maskW(:,:,k-1) + sqrt(Nsqrt(:,:,k))*dzt(k)*maskT(:,:,k) 
  enddo


  ! calculate phi_n    \sim cos(xi/c_n)*N^0.5
  phin(:,:,:,1) = 1.0
  phinW(:,:,:,1) = 1.0
  do n=2,nmodes
   do k=1,nz   
    phin(:,:,k,n)  = cos(   xi(:,:,k)/max(1d-22,cn(:,:,n)) )*Nsqrt(:,:,k)**0.25*maskT(:,:,k)
    phinW(:,:,k,n) = cos(  xiW(:,:,k)/max(1d-22,cnW(:,:,n)) )*Nsqr_lim(:,:,k)**0.25 *maskW(:,:,k)
   enddo
  enddo




  ! calculate (f/N)dphi_n/dz  \sim sin(xi/c_n)*N^0.5
  phizn(:,:,:,1)  = 0.0
  phiznW(:,:,:,1) = 0.0
  do n=2,nmodes
   do k=1,nz   
    phizn(:,:,k,n)  = sin(   xi(:,:,k)/max(1d-22,cn(:,:,n)) )*Nsqrt(:,:,k)**0.25*maskT(:,:,k)
    phiznW(:,:,k,n) = sin(  xiW(:,:,k)/max(1d-22,cnW(:,:,n)) )*Nsqr_lim(:,:,k)**0.25 *maskW(:,:,k)
   enddo
  enddo


  ! normalisation of phi_n  and  (f/N)dphi_n/dz 
  do n=1,nmodes
   ! int_(-h)^0 dz phi_n^2  = 1 
   norm=0
   do k=1,nz   
    norm = norm+ phin(:,:,k,n)**2*dzt(k)*maskT(:,:,k)
   enddo
   do k=1,nz   
     where( norm/=0) phin(:,:,k,n) = phin(:,:,k,n)/norm**0.5
   enddo
   ! int_(-h)^0  ( (f/N)dphi_n/dz)^2 dz  = 1/R^2_n
   norm=0
   do k=1,nz   
    norm = norm+ phizn(:,:,k,n)**2*dzt(k)*maskT(:,:,k)
   enddo
   do k=1,nz   
     where( norm/=0) phizn(:,:,k,n) = phizn(:,:,k,n)/norm**0.5/max(1d-12,Rn(:,:,n))
   enddo
  enddo

  ! also on W grid
  do n=1,nmodes
   norm=0
   do k=1,nz   
    norm = norm+ phinW(:,:,k,n)**2*dzw(k)*maskW(:,:,k)
   enddo
   do k=1,nz   
     where( norm/=0) phinW(:,:,k,n) = phinW(:,:,k,n)/norm**0.5
   enddo
   norm=0
   do k=1,nz   
    norm = norm+ phiznW(:,:,k,n)**2*dzw(k)*maskW(:,:,k)
   enddo
   do k=1,nz   
     where( norm/=0) phiznW(:,:,k,n) = phiznW(:,:,k,n)/norm**0.5/max(1d-12,Rn(:,:,n))
   enddo
  enddo


end subroutine calc_rossmix_vertical_struct_fct




