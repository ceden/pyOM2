



subroutine biharmonic_thickness_explicit
 use main_module
 use biharmonic_thickness_module
 ! p_t u = ... - (K (nabla^2 u)_z)_z, solve with subcycling
 implicit none
 integer :: i,j,k,n
 real*8 :: del2(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz),fxa,fxb,fxc,f_loc,Nsqr_loc!,f2
 real*8 :: aloc(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz),dt_loc
 real*8 :: kappa(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) 
 real*8 :: visc(js_pe-onx:je_pe+onx)
 !real*8, dimension(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz)  :: c1,c2,c5,c6
 !real*8, dimension(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx)  :: c3,c4
 real*8 :: c5(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz)
 real*8, save, allocatable, dimension(:,:,:)  :: c1_u,c2_u,c6_u
 real*8, save, allocatable, dimension(:,:)    :: c3_u,c4_u
 real*8, save, allocatable, dimension(:,:,:)  :: c1_v,c2_v,c6_v
 real*8, save, allocatable, dimension(:,:)    :: c3_v,c4_v
 
 logical, save :: first = .true. 
 

 dt_loc = dt_mom/biharmonic_thickness_mixing_iter


 if (first) then
  allocate( c1_u(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) )
  allocate( c2_u(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) )
  allocate( c6_u(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) )
  allocate( c3_u(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx) )
  allocate( c4_u(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx) )
  allocate( c1_v(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) )
  allocate( c2_v(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) )
  allocate( c6_v(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) )
  allocate( c3_v(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx) )
  allocate( c4_v(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx) )
 endif
 
 !---------------------------------------------------------------------------------
 ! vertical viscosity for thickness mixing
 !---------------------------------------------------------------------------------
 visc = A_thkbi * cost**(2*biha_friction_cosPower)
 do k=1,nz
  do j=js_pe-onx,je_pe+onx  
   fxc = A_thkbi_cut
   if (enable_biharmonic_thickness_scale_A_thkbi_cut)  then
      fxc = fxc*(tanh((yt(j)+60)/5.)+1)/2.*(1-tanh((yt(j)-60)/5.))/2.
     !fxc = fxc*cost(j)**(2*biha_friction_cosPower)
   endif 
   do i=is_pe-onx,ie_pe+onx    
    f_loc = max(thkbi_f_min,abs(coriolis_t(i,j)) )
    Nsqr_loc = max(1d-24,Nsqr(i,j,k,tau))
    
    fxb = visc(j) * f_loc*f_loc/Nsqr_loc * maskW(i,j,k)
    !fxa = (tanh(( sqrt(Nsqr_loc)/f_loc - 10.)*2)+1.)/2. ! clipping function for small N
    ! upper threshold from stability criterium
    kappa(i,j,k) = min(  fxb, fxc*max(cost(j)*dxt(i),dyt(j))**2*dzt(k)**2/(16*dt_loc) ) 
    !kappa(i,j,k) = fxb 
    
    !if (  fxb/max(cost(j)*dxt(i),dyt(j))**2/dzt(k)**2*16*dt_loc>1 ) then
    ! print*,' EOF violated at ',i,j,k
    ! call halt_stop(' ')
    !endif
    K_thkbi(i,j,k) = kappa(i,j,k)/max(1e-12,fxb)

   enddo
  enddo
 enddo

 if (enable_biharmonic_thickness_cut_upper) then
   do k=1,nz
    if (abs(zt(k))<A_thkbi_cut_upper_len) kappa(:,:,k) = 0d0
   enddo
 endif
 
 !---------------------------------------------------------------------------------
 ! Zonal velocity
 !---------------------------------------------------------------------------------

 aloc = u(:,:,:,tau)
 
 if (first) then
 ! precalculate coefficients to speed up 
  do j=js_pe-2,je_pe+1   
   do i=is_pe-2,ie_pe+1
    c1_u(i,j,:) = 1./(cost(j)*dxt(i+1))*maskU(i+1,j,:)*maskU(i,j,:)
    c2_u(i,j,:) = 1./dyu(j)*maskU(i,j+1,:)*maskU(i,j,:)*cosu(j)
    c3_u(i,j) = 1./(cost(j)*dxu(i)) 
    c4_u(i,j) = 1./(cost(j)*dyt(j)) 
   enddo
  enddo
  do k=1,nz
   c6_u(:,:,k) = dt_loc/dzt(k)*maskU(:,:,k)
  enddo
 endif
 
 
 do k=1,nz-1
   do j=js_pe-1,je_pe
    do i=is_pe-1,ie_pe 
     fxa = 0.5*(kappa(i,j,k)+kappa(i+1,j,k))
     c5(i,j,k) = fxa/dzw(k)*maskU(i,j,k+1)*maskU(i,j,k)
    enddo
   enddo
 enddo

 do n=1,biharmonic_thickness_mixing_iter ! sub cycling within time step
 
  call border_exchg_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,aloc) 
  call setcyclic_xyz   (is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,aloc)
  !call set_obc_boundary_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,aloc) 

  do j=js_pe-1,je_pe+1   
   do i=is_pe-2,ie_pe+1
    flux_east(i,j,:) = (aloc(i+1,j,:)-aloc(i,j,:))*c1_u(i,j,:) ! /(cost(j)*dxt(i+1))*maskU(i+1,j,:)*maskU(i,j,:)
   enddo
  enddo
  do j=js_pe-2,je_pe+1
   do i=is_pe-1,ie_pe+1
    flux_north(i,j,:) = (aloc(i,j+1,:)-aloc(i,j,:))*c2_u(i,j,:) !/dyu(j)*maskU(i,j+1,:)*maskU(i,j,:)*cosu(j)
   enddo 
  enddo 

  do j=js_pe-1,je_pe+1
   do i=is_pe-1,ie_pe+1
    del2(i,j,:) = (flux_east(i,j,:) - flux_east(i-1,j,:))*c3_u(i,j) & !/(cost(j)*dxu(i)) &
                 +(flux_north(i,j,:) - flux_north(i,j-1,:))*c4_u(i,j)  !/(cost(j)*dyt(j)) 
   enddo
  enddo

  do k=1,nz-1
   do j=js_pe-1,je_pe
    do i=is_pe-1,ie_pe 
     !fxa = 0.5*(kappa(i,j,k)+kappa(i+1,j,k))
     !flux_top(i,j,k) = fxa*(del2(i,j,k+1)-del2(i,j,k))/dzw(k)*maskU(i,j,k+1)*maskU(i,j,k)
     flux_top(i,j,k) = (del2(i,j,k+1)-del2(i,j,k))*c5(i,j,k)
    enddo
   enddo
  enddo
  flux_top(:,:,nz)=0d0
  
  !k=1; aloc(:,:,k) = aloc(:,:,k) - dt_loc*flux_top(:,:,k)/dzt(k)*maskU(:,:,k)
  k=1; aloc(:,:,k) = aloc(:,:,k) - flux_top(:,:,k)*c6_u(:,:,k)
  do k=2,nz
   !aloc(:,:,k) = aloc(:,:,k) - dt_loc*(flux_top(:,:,k)-flux_top(:,:,k-1))/dzt(k)*maskU(:,:,k)
   aloc(:,:,k) = aloc(:,:,k) - (flux_top(:,:,k)-flux_top(:,:,k-1))*c6_u(:,:,k)
  enddo
  
 enddo
 
 aloc = (aloc-u(:,:,:,tau))/dt_mom
 du_mix = du_mix + aloc
  
 if (enable_conserve_energy) then
  !---------------------------------------------------------------------------------
  ! dissipation
  !---------------------------------------------------------------------------------  
  K_diss_gm = 0d0
  aloc = -maskU*u(:,:,:,tau)*aloc
  call calc_diss(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,aloc,K_diss_gm,'U')
 endif


 !---------------------------------------------------------------------------------
 ! Meridional velocity
 !---------------------------------------------------------------------------------

 aloc = v(:,:,:,tau)
 
 ! precalculate coefficients to speed up 
 if (first) then
  do j=js_pe-2,je_pe+1
   do i=is_pe-2,ie_pe+1
    c1_v(i,j,:) = 1./(cosu(j)*dxu(i)) *maskV(i+1,j,:)*maskV(i,j,:)
    c2_v(i,j,:) = 1./dyt(j+1)*cost(j+1)*maskV(i,j,:)*maskV(i,j+1,:)
    c3_v(i,j) =  1./(cosu(j)*dxt(i)) 
    c4_v(i,j) =  1./(dyu(j)*cosu(j)) 
   enddo
  enddo
  do k=1,nz
   c6_v(:,:,k) = dt_loc/dzt(k)*maskV(:,:,k)
  enddo
 endif
 
 do k=1,nz-1
   do j=js_pe-1,je_pe
    do i=is_pe-1,ie_pe 
     fxa = 0.5*(kappa(i,j,k)+kappa(i,j+1,k))          
     c5(i,j,k) = fxa/dzw(k)*maskV(i,j,k+1)*maskV(i,j,k)
    enddo
   enddo
 enddo 

 do n=1,biharmonic_thickness_mixing_iter ! sub cycling within time step

  call border_exchg_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,aloc) 
  call setcyclic_xyz   (is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,aloc)
  !call set_obc_boundary_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,aloc) 

  do j=js_pe-1,je_pe+1
   do i=is_pe-2,ie_pe+1
    flux_east(i,j,:) = (aloc(i+1,j,:)-aloc(i,j,:))*c1_v(i,j,:)!/(cosu(j)*dxu(i)) *maskV(i+1,j,:)*maskV(i,j,:)
   enddo
  enddo
  do j=js_pe-2,je_pe+1
   do i=is_pe-1,ie_pe+1
    flux_north(i,j,:) = (aloc(i,j+1,:)-aloc(i,j,:) )*c2_v(i,j,:)!/dyt(j+1)*cost(j+1)*maskV(i,j,:)*maskV(i,j+1,:)
   enddo
  enddo

  do j=js_pe-1,je_pe+1
   do i=is_pe-1,ie_pe+1
    del2(i,j,:) =  (flux_east(i,j,:) - flux_east(i-1,j,:))*c3_v(i,j) & !/(cosu(j)*dxt(i))  &
                  +(flux_north(i,j,:) - flux_north(i,j-1,:))*c4_v(i,j)!/(dyu(j)*cosu(j)) 
   enddo
  enddo


  do k=1,nz-1
   do j=js_pe-1,je_pe
    do i=is_pe-1,ie_pe 
     !fxa = 0.5*(kappa(i,j,k)+kappa(i,j+1,k))          
     !flux_top(i,j,k) = fxa*(del2(i,j,k+1)-del2(i,j,k))/dzw(k)*maskV(i,j,k+1)*maskV(i,j,k)
     flux_top(i,j,k) = (del2(i,j,k+1)-del2(i,j,k))*c5(i,j,k)
    enddo
   enddo
  enddo
  flux_top(:,:,nz)=0d0
  
  !k=1; aloc(:,:,k) = aloc(:,:,k) - dt_loc*flux_top(:,:,k)/dzt(k)*maskV(:,:,k)
  k=1; aloc(:,:,k) = aloc(:,:,k) - flux_top(:,:,k)*c6_v(:,:,k)
  do k=2,nz
   !aloc(:,:,k) = aloc(:,:,k) - dt_loc*(flux_top(:,:,k)-flux_top(:,:,k-1))/dzt(k)*maskV(:,:,k)
   aloc(:,:,k) = aloc(:,:,k) - (flux_top(:,:,k)-flux_top(:,:,k-1))*c6_v(:,:,k)
  enddo

 enddo
 
 aloc = (aloc-v(:,:,:,tau))/dt_mom
 dv_mix = dv_mix + aloc

 if (enable_conserve_energy) then
  !---------------------------------------------------------------------------------
  ! dissipation
  !---------------------------------------------------------------------------------
  aloc = -maskV*v(:,:,:,tau)*aloc
  call calc_diss(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,aloc,K_diss_gm,'V')
 endif

 first = .false.
end subroutine biharmonic_thickness_explicit





