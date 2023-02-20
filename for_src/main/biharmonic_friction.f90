
  

subroutine biharmonic_friction
!=======================================================================
! horizontal biharmonic friction   
! dissipation is calculated and added to K_diss_h
!=======================================================================
 use main_module   
 implicit none
 integer :: i,j
 real*8 :: del2(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz),fxa
 real*8 :: diss(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz),visc

 if (.not.enable_hydrostatic) call halt_stop('biharmonic mixing for non-hydrostatic not yet implemented')

 visc = sqrt(abs(A_hbi))
 !---------------------------------------------------------------------------------
 ! Zonal velocity
 !---------------------------------------------------------------------------------
 do j=js_pe-1,je_pe+1
   fxa = cost(j)**biha_friction_cosPower
   do i=is_pe-2,ie_pe+1
    flux_east(i,j,:)=visc*fxa*(u(i+1,j,:,tau)-u(i,j,:,tau))/(cost(j)*dxt(i+1))*maskU(i+1,j,:)*maskU(i,j,:)
   enddo
 enddo
 do j=js_pe-2,je_pe+1
  fxa = cosu(j)**biha_friction_cosPower
  do i=is_pe-1,ie_pe+1
   flux_north(i,j,:)=visc*fxa*(u(i,j+1,:,tau)-u(i,j,:,tau))/dyu(j)*maskU(i,j+1,:)*maskU(i,j,:)*cosu(j)
  enddo 
 enddo 

 do j=js_pe-1,je_pe+1
   do i=is_pe-1,ie_pe+1
    del2(i,j,:)= (flux_east(i,j,:) - flux_east(i-1,j,:))/(cost(j)*dxu(i)) &
                +(flux_north(i,j,:) - flux_north(i,j-1,:))/(cost(j)*dyt(j)) 
  enddo
 enddo

 do j=js_pe,je_pe
   fxa = cost(j)**biha_friction_cosPower
   do i=is_pe-1,ie_pe
    flux_east(i,j,:)=visc*fxa*(del2(i+1,j,:)-del2(i,j,:))/(cost(j)*dxt(i+1))*maskU(i+1,j,:)*maskU(i,j,:)
   enddo
 enddo
 do j=js_pe-1,je_pe
   fxa = cosu(j)**biha_friction_cosPower
   do i=is_pe,ie_pe
    flux_north(i,j,:)=visc*fxa*(del2(i,j+1,:)-del2(i,j,:))/dyu(j)*maskU(i,j+1,:)*maskU(i,j,:)*cosu(j)
   enddo 
 enddo 

 !---------------------------------------------------------------------------------
 ! update tendency 
 !---------------------------------------------------------------------------------
 do j=js_pe,je_pe
   do i=is_pe,ie_pe
    du_mix(i,j,:)= du_mix(i,j,:) - maskU(i,j,:)*((flux_east(i,j,:) - flux_east(i-1,j,:))/(cost(j)*dxu(i)) &
                                                +(flux_north(i,j,:) - flux_north(i,j-1,:))/(cost(j)*dyt(j)) )
  enddo
 enddo

 if (enable_conserve_energy) then
 !---------------------------------------------------------------------------------
 ! diagnose dissipation by biharmonic friction of u
 !---------------------------------------------------------------------------------
  do j=js_pe,je_pe
   fxa = cost(j)**biha_friction_cosPower
   do i=is_pe-2,ie_pe
    flux_east(i,j,:)=visc*fxa*(u(i+1,j,:,tau)-u(i,j,:,tau))/(cost(j)*dxt(i+1))*maskU(i+1,j,:)*maskU(i,j,:)
   enddo
  enddo
  do j=js_pe-1,je_pe
   fxa = cosu(j)**biha_friction_cosPower
   do i=is_pe-1,ie_pe
    flux_north(i,j,:)=visc*fxa*(u(i,j+1,:,tau)-u(i,j,:,tau))/dyu(j)*maskU(i,j+1,:)*maskU(i,j,:)*cosu(j)
   enddo 
  enddo 
  do j=js_pe,je_pe
   do i=is_pe-1,ie_pe
    diss(i,j,:)=  maskU(i,j,:)*((flux_east(i,j,:) - flux_east(i-1,j,:))/(cost(j)*dxu(i)) &
                              +(flux_north(i,j,:) - flux_north(i,j-1,:))/(cost(j)*dyt(j)) )**2
   enddo
  enddo
  call calc_diss(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,diss,K_diss_h,'U')
 endif

 !---------------------------------------------------------------------------------
 ! Meridional velocity
 !---------------------------------------------------------------------------------
 do j=js_pe-1,je_pe+1
   fxa = cosu(j)**biha_friction_cosPower
   do i=is_pe-2,ie_pe+1
    flux_east(i,j,:)=visc*fxa*(v(i+1,j,:,tau)-v(i,j,:,tau))/(cosu(j)*dxu(i)) *maskV(i+1,j,:)*maskV(i,j,:)
   enddo
 enddo
 do j=js_pe-2,je_pe+1
   fxa = cost(j+1)**biha_friction_cosPower
   do i=is_pe-1,ie_pe+1
    flux_north(i,j,:)=visc*fxa*(v(i,j+1,:,tau)-v(i,j,:,tau) )/dyt(j+1)*cost(j+1)*maskV(i,j,:)*maskV(i,j+1,:)
   enddo
 enddo

 do j=js_pe-1,je_pe+1
  do i=is_pe-1,ie_pe+1
    del2(i,j,:)=  (flux_east(i,j,:) - flux_east(i-1,j,:))/(cosu(j)*dxt(i))  &
                 +(flux_north(i,j,:) - flux_north(i,j-1,:))/(dyu(j)*cosu(j)) 
  enddo
 enddo

 do j=js_pe,je_pe
   fxa = cosu(j)**biha_friction_cosPower
   do i=is_pe-1,ie_pe
    flux_east(i,j,:)=visc*fxa*(del2(i+1,j,:)-del2(i,j,:))/(cosu(j)*dxu(i)) *maskV(i+1,j,:)*maskV(i,j,:)
   enddo
 enddo
 do j=js_pe-1,je_pe
   fxa = cost(j+1)**biha_friction_cosPower
   do i=is_pe,ie_pe
    flux_north(i,j,:)=visc*fxa*(del2(i,j+1,:)-del2(i,j,:) )/dyt(j+1)*cost(j+1)*maskV(i,j,:)*maskV(i,j+1,:)
   enddo
 enddo

 !---------------------------------------------------------------------------------
 ! update tendency 
 !---------------------------------------------------------------------------------
 do j=js_pe,je_pe
   do i=is_pe,ie_pe
    dv_mix(i,j,:)= dv_mix(i,j,:) - maskV(i,j,:)*( (flux_east(i,j,:) - flux_east(i-1,j,:))/(cosu(j)*dxt(i))  &
                                                 +(flux_north(i,j,:) - flux_north(i,j-1,:))/(dyu(j)*cosu(j)) )
  enddo
 enddo

 if (enable_conserve_energy) then
 !---------------------------------------------------------------------------------
 ! diagnose dissipation by biharmonic friction of v
 !---------------------------------------------------------------------------------
  do j=js_pe-1,je_pe
   fxa = cosu(j)**biha_friction_cosPower
   do i=is_pe-1,ie_pe
    flux_east(i,j,:)=visc*fxa*(v(i+1,j,:,tau)-v(i,j,:,tau))/(cosu(j)*dxu(i)) *maskV(i+1,j,:)*maskV(i,j,:)
   enddo
  enddo
  do j=js_pe-2,je_pe
   fxa = cost(j+1)**biha_friction_cosPower
   do i=is_pe,ie_pe
    flux_north(i,j,:)=visc*fxa*(v(i,j+1,:,tau)-v(i,j,:,tau) )/dyt(j+1)*cost(j+1)*maskV(i,j,:)*maskV(i,j+1,:)
   enddo
  enddo
  do j=js_pe-1,je_pe
   do i=is_pe,ie_pe
    diss(i,j,:)= maskV(i,j,:)*( (flux_east(i,j,:) - flux_east(i-1,j,:))/(cosu(j)*dxt(i))  &
                               +(flux_north(i,j,:) - flux_north(i,j-1,:))/(dyu(j)*cosu(j)) )**2
   enddo
  enddo    
  call calc_diss(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,diss,K_diss_h,'V')
 endif

end subroutine biharmonic_friction





subroutine biharmonic_thickness_mixing
 use main_module
 use isoneutral_module 
 implicit none
 integer :: i,j,k
 real*8 :: del2(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz),fxa,fxb,f_loc,Nsqr_loc!,f2
 real*8 :: diss(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz)
 real*8 :: kappa(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz)

 if (.not.enable_TEM_friction .and. enable_conserve_energy) K_diss_gm = 0d0

 do k=1,nz
  do j=js_pe-onx,je_pe+onx   
   do i=is_pe-onx,ie_pe+onx    
    f_loc = max(thkbi_f_min,abs(coriolis_t(i,j)) )
    Nsqr_loc = max(1d-12,Nsqr(i,j,k,tau))
    fxa = (tanh(( sqrt(Nsqr_loc)/f_loc - 10.)*2)+1.)/2.
    fxb = A_thkbi * (cost(j)**biha_friction_cosPower)**2
    kappa(i,j,k) = fxa*fxb*f_loc**2/Nsqr_loc*maskT(i,j,k)
    ! different clipping function
    !f1 = max(thkbi_f_min, abs(coriolis_t(i  ,j)))
    !f2 = max(thkbi_f_min, abs(coriolis_t(i+1,j)))
    !kappa(i,j,k) =fxb*(  0.5*min(0.01d0,f1**2/max(1d-12,Nsqr(i  ,j,k,tau)) ) &
    !                   + 0.5*min(0.01d0,f2**2/max(1d-12,Nsqr(i+1,j,k,tau)) ) )        
   enddo
  enddo
 enddo

 !---------------------------------------------------------------------------------
 ! Zonal velocity
 !---------------------------------------------------------------------------------
 do j=js_pe-1,je_pe+1   
   do i=is_pe-2,ie_pe+1
    flux_east(i,j,:)=(u(i+1,j,:,tau)-u(i,j,:,tau))/(cost(j)*dxt(i+1))*maskU(i+1,j,:)*maskU(i,j,:)
   enddo
 enddo
 do j=js_pe-2,je_pe+1
  do i=is_pe-1,ie_pe+1
   flux_north(i,j,:)=(u(i,j+1,:,tau)-u(i,j,:,tau))/dyu(j)*maskU(i,j+1,:)*maskU(i,j,:)*cosu(j)
  enddo 
 enddo 

 do j=js_pe-1,je_pe+1
   do i=is_pe-1,ie_pe+1
    del2(i,j,:)= (flux_east(i,j,:) - flux_east(i-1,j,:))/(cost(j)*dxu(i)) &
                +(flux_north(i,j,:) - flux_north(i,j-1,:))/(cost(j)*dyt(j)) 
  enddo
 enddo

 do k=1,nz-1
  do j=js_pe-1,je_pe
   do i=is_pe-1,ie_pe 
    fxa = 0.5*(kappa(i,j,k)+kappa(i+1,j,k))
    flux_top(i,j,k)=fxa*(del2(i,j,k+1)-del2(i,j,k))/dzw(k)*maskU(i,j,k+1)*maskU(i,j,k)
   enddo
  enddo
 enddo
 flux_top(:,:,nz)=0d0
 k=1; du_mix(:,:,k) = du_mix(:,:,k) - flux_top(:,:,k)/dzt(k)*maskU(:,:,k)
 do k=2,nz
   du_mix(:,:,k) = du_mix(:,:,k) - (flux_top(:,:,k)-flux_top(:,:,k-1))/dzt(k)*maskU(:,:,k)
 enddo
 
 if (enable_conserve_energy) then
 !---------------------------------------------------------------------------------
 ! dissipation
 !---------------------------------------------------------------------------------
  k=1; diss(:,:,k) = maskU(:,:,k)*u(:,:,k,tau)*flux_top(:,:,k)/dzt(k)*maskU(:,:,k)
  do k=2,nz
   diss(:,:,k) = maskU(:,:,k)*u(:,:,k,tau)*(flux_top(:,:,k)-flux_top(:,:,k-1))/dzt(k)*maskU(:,:,k)
  enddo  
  call calc_diss(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,diss,K_diss_gm,'U')
 endif


 !---------------------------------------------------------------------------------
 ! Meridional velocity
 !---------------------------------------------------------------------------------
 do j=js_pe-1,je_pe+1
   do i=is_pe-2,ie_pe+1
    flux_east(i,j,:)=(v(i+1,j,:,tau)-v(i,j,:,tau))/(cosu(j)*dxu(i)) *maskV(i+1,j,:)*maskV(i,j,:)
   enddo
 enddo
 do j=js_pe-2,je_pe+1
   do i=is_pe-1,ie_pe+1
    flux_north(i,j,:)=(v(i,j+1,:,tau)-v(i,j,:,tau) )/dyt(j+1)*cost(j+1)*maskV(i,j,:)*maskV(i,j+1,:)
   enddo
 enddo

 do j=js_pe-1,je_pe+1
  do i=is_pe-1,ie_pe+1
    del2(i,j,:)=  (flux_east(i,j,:) - flux_east(i-1,j,:))/(cosu(j)*dxt(i))  &
                 +(flux_north(i,j,:) - flux_north(i,j-1,:))/(dyu(j)*cosu(j)) 
  enddo
 enddo

 do k=1,nz-1
  do j=js_pe-1,je_pe
   do i=is_pe-1,ie_pe 
    fxa = 0.5*(kappa(i,j,k)+kappa(i,j+1,k))          
    flux_top(i,j,k)=fxa*(del2(i,j,k+1)-del2(i,j,k))/dzw(k)*maskV(i,j,k+1)*maskV(i,j,k)
   enddo
  enddo
 enddo
 flux_top(:,:,nz)=0d0
 k=1; dv_mix(:,:,k) = dv_mix(:,:,k) - flux_top(:,:,k)/dzt(k)*maskV(:,:,k)
 do k=2,nz
   dv_mix(:,:,k) = dv_mix(:,:,k) - (flux_top(:,:,k)-flux_top(:,:,k-1))/dzt(k)*maskV(:,:,k)
 enddo

 if (enable_conserve_energy) then
 !---------------------------------------------------------------------------------
 ! dissipation
 !---------------------------------------------------------------------------------
  k=1; diss(:,:,k) = maskV(:,:,k)*v(:,:,k,tau)*flux_top(:,:,k)/dzt(k)*maskV(:,:,k)
  do k=2,nz
   diss(:,:,k) = maskV(:,:,k)*v(:,:,k,tau)*(flux_top(:,:,k)-flux_top(:,:,k-1))/dzt(k)*maskV(:,:,k)
  enddo
  call calc_diss(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,diss,K_diss_gm,'V')
 endif

end subroutine biharmonic_thickness_mixing




