 



subroutine advect_tracer(is_,ie_,js_,je_,nz_,tr,dtr)
!=======================================================================
! calculate time tendency of a tracer due to advection
!=======================================================================
 use main_module   
 implicit none
 integer, intent(in) :: is_,ie_,js_,je_,nz_
 real*8, intent(inout) :: dtr(is_:ie_,js_:je_,nz_),tr(is_:ie_,js_:je_,nz_)
 integer :: i,j,k
 if (enable_superbee_advection) then
  call adv_flux_superbee(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,flux_east,flux_north,flux_top,tr)
 elseif (enable_upwind3_advection) then
  call adv_flux_upwind3(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,flux_east,flux_north,flux_top,tr)
 elseif (enable_dst3_advection) then
  call adv_flux_dst3(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,flux_east,flux_north,flux_top,tr)
 else
  call adv_flux_2nd(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,flux_east,flux_north,flux_top,tr)
 endif
 do j=js_pe,je_pe
   do i=is_pe,ie_pe
      dtr(i,j,:)=maskT(i,j,:)* (-( flux_east(i,j,:)-  flux_east(i-1,j,:))/(cost(j)*dxt(i)) &
                                -(flux_north(i,j,:)- flux_north(i,j-1,:))/(cost(j)*dyt(j)) )
   enddo
 enddo
 k=1; dtr(:,:,k)=dtr(:,:,k)-maskT(:,:,k)*flux_top(:,:,k)/dzt(k)   
 do k=2,nz
   dtr(:,:,k)=dtr(:,:,k)-maskT(:,:,k)*(flux_top(:,:,k)- flux_top(:,:,k-1))/dzt(k)   
 enddo
end subroutine advect_tracer




!subroutine advect_temperature
!=======================================================================
! integrate temperature 
!=======================================================================
! use main_module   
! implicit none
! integer :: i,j,k
! 
! if (enable_superbee_advection) then
!  call adv_flux_superbee(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,flux_east,flux_north,flux_top,temp(:,:,:,tau))
! else
!  call adv_flux_2nd(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,flux_east,flux_north,flux_top,temp(:,:,:,tau))
! endif
! do j=js_pe,je_pe
!   do i=is_pe,ie_pe
!      dtemp(i,j,:,tau)=maskT(i,j,:)* (-( flux_east(i,j,:)-  flux_east(i-1,j,:))/(cost(j)*dxt(i)) &
!                                      -(flux_north(i,j,:)- flux_north(i,j-1,:))/(cost(j)*dyt(j)) )
!   enddo
! enddo
! k=1; dtemp(:,:,k,tau)=dtemp(:,:,k,tau)-maskT(:,:,k)*flux_top(:,:,k)/dzt(k)   
! do k=2,nz
!   dtemp(:,:,k,tau)=dtemp(:,:,k,tau)-maskT(:,:,k)*(flux_top(:,:,k)- flux_top(:,:,k-1))/dzt(k)   
! enddo
!end subroutine advect_temperature



!subroutine advect_salinity
!=======================================================================
! integrate salinity
!=======================================================================
! use main_module   
! implicit none
! integer :: i,j,k
! if (enable_superbee_advection) then
!  call adv_flux_superbee(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,flux_east,flux_north,flux_top,salt(:,:,:,tau))
! else
!  call adv_flux_2nd(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,flux_east,flux_north,flux_top,salt(:,:,:,tau))
! endif
! do j=js_pe,je_pe
!   do i=is_pe,ie_pe
!      dsalt(i,j,:,tau)=maskT(i,j,:)* (-( flux_east(i,j,:)-  flux_east(i-1,j,:))/(cost(j)*dxt(i)) &
!                                      -(flux_north(i,j,:)- flux_north(i,j-1,:))/(cost(j)*dyt(j)) )
!   enddo
! enddo
! k=1; dsalt(:,:,k,tau)=dsalt(:,:,k,tau)-maskT(:,:,k)*flux_top(:,:,k)/dzt(k)   
! do k=2,nz
!   dsalt(:,:,k,tau)=dsalt(:,:,k,tau)-maskT(:,:,k)*(flux_top(:,:,k)- flux_top(:,:,k-1))/dzt(k)   
! enddo
!end subroutine advect_salinity






subroutine calculate_velocity_on_wgrid
!---------------------------------------------------------------------------------
! calculates advection velocity for tracer on W grid
!---------------------------------------------------------------------------------
 use main_module
 implicit none
 integer :: i,j,k
 !real*8 :: fxa,fxb

 ! lateral advection velocities on W grid
 do k=1,nz-1
  u_wgrid(:,:,k) = u(:,:,k+1,tau)*maskU(:,:,k+1)*0.5*dzt(k+1)/dzw(k) + u(:,:,k,tau)*maskU(:,:,k)*0.5*dzt(k)/dzw(k)
  v_wgrid(:,:,k) = v(:,:,k+1,tau)*maskV(:,:,k+1)*0.5*dzt(k+1)/dzw(k) + v(:,:,k,tau)*maskV(:,:,k)*0.5*dzt(k)/dzw(k)
 enddo
 k=nz
 u_wgrid(:,:,k) = u(:,:,k,tau)*maskU(:,:,k)*0.5*dzt(k)/dzw(k)
 v_wgrid(:,:,k) = v(:,:,k,tau)*maskV(:,:,k)*0.5*dzt(k)/dzw(k)

 ! redirect velocity at bottom and at topography
 k=1
 u_wgrid(:,:,k) = u_wgrid(:,:,k) + u(:,:,k,tau)*maskU(:,:,k)*0.5*dzt(k)/dzw(k)
 v_wgrid(:,:,k) = v_wgrid(:,:,k) + v(:,:,k,tau)*maskV(:,:,k)*0.5*dzt(k)/dzw(k) 
 do k=1,nz-1
  do j=js_pe-onx,je_pe+onx
    do i=is_pe-onx,ie_pe+onx-1
      if (maskW(i,j,k)*maskW(i+1,j,k) == 0d0 ) then
         u_wgrid(i,j,k+1) = u_wgrid(i,j,k+1)+u_wgrid(i,j,k)*dzw(k)/dzw(k+1)
         u_wgrid(i,j,k) = 0d0
      endif
   enddo
  enddo
  do j=js_pe-onx,je_pe+onx-1
    do i=is_pe-onx,ie_pe+onx
      if (maskW(i,j,k)*maskW(i,j+1,k) == 0d0 ) then
         v_wgrid(i,j,k+1) = v_wgrid(i,j,k+1)+v_wgrid(i,j,k)*dzw(k)/dzw(k+1)
         v_wgrid(i,j,k) = 0d0
      endif
   enddo
  enddo
 enddo

 ! vertical advection velocity on W grid from continuity
 w_wgrid(:,:,1)=0d0
 do k=1,nz
  do j=js_pe-onx+1,je_pe+onx
    do i=is_pe-onx+1,ie_pe+onx
      w_wgrid(i,j,k) = w_wgrid(i,j,max(1,k-1))-dzw(k)* &
               ((        u_wgrid(i,j,k)          -u_wgrid(i-1,j,k))/(cost(j)*dxt(i)) &
               +(cosu(j)*v_wgrid(i,j,k)-cosu(j-1)*v_wgrid(i,j-1,k))/(cost(j)*dyt(j)) )
   enddo
  enddo
 enddo

 ! test continuity
 !if ( modulo(itt*dt_tracer,ts_monint) < dt_tracer .and. .false.) then
 ! fxa=0;fxb=0;
 ! do j=js_pe,je_pe
 !  do i=is_pe,ie_pe
 !    fxa = fxa + w_wgrid(i,j,nz) *area_t(i,j)
 !    fxb = fxb +   w(i,j,nz,tau) *area_t(i,j)
 !  enddo
 ! enddo
 ! call global_sum(fxa); call global_sum(fxb); 
 ! if (my_pe==0) print'(a,e12.6,a)',' transport at sea surface on t grid = ',fxb,' m^3/s'
 ! if (my_pe==0) print'(a,e12.6,a)',' transport at sea surface on w grid = ',fxa,' m^3/s'
!
!
!  fxa=0;fxb=0;
!  do j=js_pe,je_pe
!   do i=is_pe,ie_pe
!     fxa = fxa + w_wgrid(i,j,nz)**2 *area_t(i,j)
!     fxb = fxb +   w(i,j,nz,tau)**2 *area_t(i,j)
!   enddo
!  enddo
!  call global_sum(fxa); call global_sum(fxb); 
!  if (my_pe==0) print'(a,e12.6,a)',' w variance on t grid = ',fxb,' (m^3/s)^2'
!  if (my_pe==0) print'(a,e12.6,a)',' w variance on w grid = ',fxa,' (m^3/s)^2'
!
! endif

end subroutine calculate_velocity_on_wgrid







