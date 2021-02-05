


subroutine harmonic_friction
!=======================================================================
! horizontal harmonic friction   
! dissipation is calculated and added to K_diss_h
!=======================================================================
 use main_module   
 implicit none
 integer :: i,j,k
 integer :: is,ie,js,je
 real*8 :: diss(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz),fxa

 is = is_pe-onx; ie = ie_pe+onx; js = js_pe-onx; je = je_pe+onx

 !---------------------------------------------------------------------------------
 ! Zonal velocity
 !---------------------------------------------------------------------------------
 if (enable_hor_friction_cos_scaling) then
  do j=js,je
   fxa = cost(j)**hor_friction_cosPower
   do i=is,ie-1
    flux_east(i,j,:)=fxa*A_h*(u(i+1,j,:,tau)-u(i,j,:,tau))/(cost(j)*dxt(i+1))*maskU(i+1,j,:)*maskU(i,j,:)
   enddo
  enddo
  do j=js,je-1
   fxa = cosu(j)**hor_friction_cosPower
   flux_north(:,j,:)=fxa*A_h*(u(:,j+1,:,tau)-u(:,j,:,tau))/dyu(j)*maskU(:,j+1,:)*maskU(:,j,:)*cosu(j)
  enddo 
 else
  do j=js,je
   do i=is,ie-1
    flux_east(i,j,:)=A_h*(u(i+1,j,:,tau)-u(i,j,:,tau))/(cost(j)*dxt(i+1))*maskU(i+1,j,:)*maskU(i,j,:)
   enddo
  enddo
  do j=js,je-1
    flux_north(:,j,:)=A_h*(u(:,j+1,:,tau)-u(:,j,:,tau))/dyu(j)*maskU(:,j+1,:)*maskU(:,j,:)*cosu(j)
  enddo 
 endif
 flux_east(ie,:,:)=0.
 flux_north(:,je,:)=0.

 !---------------------------------------------------------------------------------
 ! update tendency 
 !---------------------------------------------------------------------------------
 do j=js_pe,je_pe
   do i=is_pe,ie_pe
    du_mix(i,j,:)= du_mix(i,j,:) + maskU(i,j,:)*((flux_east(i,j,:) - flux_east(i-1,j,:))/(cost(j)*dxu(i)) &
                                                +(flux_north(i,j,:) - flux_north(i,j-1,:))/(cost(j)*dyt(j)) )
  enddo
 enddo

 if (enable_conserve_energy) then
 !---------------------------------------------------------------------------------
 ! diagnose dissipation by lateral friction
 !---------------------------------------------------------------------------------
  do k=1,nz
   do j=js_pe,je_pe
    do i=is_pe-1,ie_pe
     diss(i,j,k) =0.5*((u(i+1,j,k,tau)-u(i,j,k,tau))*flux_east(i,j,k) &
                      +(u(i,j,k,tau)-u(i-1,j,k,tau))*flux_east(i-1,j,k))/(cost(j)*dxu(i))  &
                 +0.5*((u(i,j+1,k,tau)-u(i,j,k,tau))*flux_north(i,j,k)+ &
                       (u(i,j,k,tau)-u(i,j-1,k,tau))*flux_north(i,j-1,k))/(cost(j)*dyt(j)) 
    enddo
   enddo
  enddo
  call calc_diss(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,diss,K_diss_h,'U')
 endif

 !---------------------------------------------------------------------------------
 ! Meridional velocity
 !---------------------------------------------------------------------------------
 if (enable_hor_friction_cos_scaling) then
  do j=js,je
   fxa = cosu(j)**hor_friction_cosPower
   do i=is,ie-1
    flux_east(i,j,:)=fxa*A_h*(v(i+1,j,:,tau)-v(i,j,:,tau))/(cosu(j)*dxu(i)) *maskV(i+1,j,:)*maskV(i,j,:)
   enddo
  enddo
  do j=js,je-1
   fxa = cost(j+1)**hor_friction_cosPower
   flux_north(:,j,:)=fxa*A_h*(v(:,j+1,:,tau)-v(:,j,:,tau) )/dyt(j+1)*cost(j+1)*maskV(:,j,:)*maskV(:,j+1,:)
  enddo
 else
  do j=js,je
   do i=is,ie-1
    flux_east(i,j,:)=A_h*(v(i+1,j,:,tau)-v(i,j,:,tau))/(cosu(j)*dxu(i)) *maskV(i+1,j,:)*maskV(i,j,:)
   enddo
  enddo
  do j=js,je-1
   flux_north(:,j,:)=A_h*(v(:,j+1,:,tau)-v(:,j,:,tau) )/dyt(j+1)*cost(j+1)*maskV(:,j,:)*maskV(:,j+1,:)
  enddo
 endif
 flux_east(ie,:,:)=0.
 flux_north(:,je,:)=0.

 !---------------------------------------------------------------------------------
 ! update tendency 
 !---------------------------------------------------------------------------------
 do j=js_pe,je_pe
   do i=is_pe,ie_pe
    dv_mix(i,j,:)= dv_mix(i,j,:) + maskV(i,j,:)*( (flux_east(i,j,:) - flux_east(i-1,j,:))/(cosu(j)*dxt(i))  &
                                                 +(flux_north(i,j,:) - flux_north(i,j-1,:))/(dyu(j)*cosu(j)) )
  enddo
 enddo

 if (enable_conserve_energy) then
 !---------------------------------------------------------------------------------
 ! diagnose dissipation by lateral friction
 !---------------------------------------------------------------------------------
  do k=1,nz
   do j=js_pe-1,je_pe
    do i=is_pe,ie_pe
     diss(i,j,k) =0.5*((v(i+1,j,k,tau)-v(i,j,k,tau))*flux_east(i,j,k)+ &
                       (v(i,j,k,tau)-v(i-1,j,k,tau))*flux_east(i-1,j,k))/(cosu(j)*dxt(i)) &
                + 0.5*((v(i,j+1,k,tau)-v(i,j,k,tau))*flux_north(i,j,k)+ &
                       (v(i,j,k,tau)-v(i,j-1,k,tau))*flux_north(i,j-1,k))/(cosu(j)*dyu(j)) 
    enddo
   enddo
  enddo
  call calc_diss(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,diss,K_diss_h,'V')
 endif

 if (.not.enable_hydrostatic) then

  if (enable_hor_friction_cos_scaling) then
   if (my_pe==0) print'(/a/)','ERROR: scaling of lateral friction for vertical velocity not implemented'
   call halt_stop(' in hamronic_friction')
  endif

  do j=js,je
   do i=is,ie-1
    flux_east(i,j,:)=A_h*(w(i+1,j,:,tau)-w(i,j,:,tau))/(cost(j)*dxu(i)) *maskW(i+1,j,:)*maskW(i,j,:)
   enddo
  enddo
  do j=js,je-1
   flux_north(:,j,:)=A_h*(w(:,j+1,:,tau)-w(:,j,:,tau))/dyu(j)*maskW(:,j+1,:)*maskW(:,j,:)*cosu(j)
  enddo
  flux_east(ie,:,:)=0.
  flux_north(:,je,:)=0.

 !---------------------------------------------------------------------------------
 ! update tendency 
 !---------------------------------------------------------------------------------
  do j=js_pe,je_pe
   do i=is_pe,ie_pe
    dw_mix(i,j,:)= dw_mix(i,j,:) + maskW(i,j,:)*( (flux_east(i,j,:) - flux_east(i-1,j,:))/(cost(j)*dxt(i))  &
                                                 +(flux_north(i,j,:) - flux_north(i,j-1,:))/(dyt(j)*cost(j)) )
  enddo
  enddo

 !---------------------------------------------------------------------------------
 ! diagnose dissipation by lateral friction
 !---------------------------------------------------------------------------------
  ! to be implemented
 endif
end subroutine harmonic_friction




