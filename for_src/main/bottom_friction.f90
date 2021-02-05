


subroutine rayleigh_friction
!=======================================================================
!  interior Rayleigh friction   
!  dissipation is calculated and added to K_diss_bot
!=======================================================================
 use main_module   
 implicit none
 integer :: k
 real*8 :: diss(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz)

 do k=1,nz
  du_mix(:,:,k)=du_mix(:,:,k) - maskU(:,:,k)*r_ray*u(:,:,k,tau)
 enddo
 if (enable_conserve_energy) then
  do k=1,nz
   diss(:,:,k) = maskU(:,:,k)*r_ray*u(:,:,k,tau)**2
  enddo
  call calc_diss(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,diss,K_diss_bot,'U')
 endif

 do k=1,nz
   dv_mix(:,:,k)=dv_mix(:,:,k) - maskV(:,:,k)*r_ray*v(:,:,k,tau)
 enddo
 if (enable_conserve_energy) then
  do k=1,nz
   diss(:,:,k) = maskV(:,:,k)*r_ray*v(:,:,k,tau)**2
  enddo
  call calc_diss(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,diss,K_diss_bot,'V')
 endif

 if (.not.enable_hydrostatic) then
  !if (my_pe==0) print'(/a/)','ERROR: rayleigh friction for vertical velocity not implemented'
  !call halt_stop(' in rayleigh_friction')
 endif
end subroutine rayleigh_friction






subroutine linear_bottom_friction
!=======================================================================
!   linear bottom friction   
!   dissipation is calculated and added to K_diss_bot
!=======================================================================
 use main_module   
 implicit none
 integer :: i,j,k
 real*8 :: diss(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz)

 if (enable_bottom_friction_var) then

 !---------------------------------------------------------------------------------
 ! with spatially varying coefficient
 !---------------------------------------------------------------------------------
  do j=js_pe,je_pe
   do i=is_pe-1,ie_pe
    k=max(kbot(i,j),kbot(i+1,j))
    if (k>0) du_mix(i,j,k)=du_mix(i,j,k) - maskU(i,j,k)*r_bot_var_u(i,j)*u(i,j,k,tau)
   enddo
  enddo
  if (enable_conserve_energy) then
   diss=0.0
   do j=js_pe,je_pe
    do i=is_pe-1,ie_pe
     k=max(kbot(i,j),kbot(i+1,j))
     if (k>0) diss(i,j,k)= maskU(i,j,k)*r_bot_var_u(i,j)*u(i,j,k,tau)**2
    enddo
   enddo
   call calc_diss(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,diss,K_diss_bot,'U')
  endif

  do j=js_pe-1,je_pe
   do i=is_pe,ie_pe
    k=max(kbot(i,j+1),kbot(i,j))
    if (k>0) dv_mix(i,j,k)=dv_mix(i,j,k) - maskV(i,j,k)*r_bot_var_v(i,j)*v(i,j,k,tau)
   enddo
  enddo
  if (enable_conserve_energy) then
   diss=0.0
   do j=js_pe-1,je_pe
    do i=is_pe,ie_pe
     k=max(kbot(i,j+1),kbot(i,j))
     if (k>0) diss(i,j,k) = maskV(i,j,k)*r_bot_var_v(i,j)*v(i,j,k,tau)**2
    enddo
   enddo
   call calc_diss(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,diss,K_diss_bot,'V')
  endif

 else
 !---------------------------------------------------------------------------------
 ! with constant coefficient
 !---------------------------------------------------------------------------------
  do j=js_pe,je_pe
   do i=is_pe-1,ie_pe
    k=max(kbot(i,j),kbot(i+1,j))
    if (k>0) du_mix(i,j,k)=du_mix(i,j,k) - maskU(i,j,k)*r_bot*u(i,j,k,tau)
   enddo
  enddo
  if (enable_conserve_energy) then
   diss=0.0
   do j=js_pe,je_pe
    do i=is_pe-1,ie_pe
     k=max(kbot(i,j),kbot(i+1,j))
     if (k>0) diss(i,j,k)= maskU(i,j,k)*r_bot*u(i,j,k,tau)**2
    enddo
   enddo
   call calc_diss(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,diss,K_diss_bot,'U')
  endif
 
  do j=js_pe-1,je_pe
   do i=is_pe,ie_pe
    k=max(kbot(i,j+1),kbot(i,j))
    if (k>0) dv_mix(i,j,k)=dv_mix(i,j,k) - maskV(i,j,k)*r_bot*v(i,j,k,tau)
   enddo
  enddo
  if (enable_conserve_energy) then
   diss=0.0
   do j=js_pe-1,je_pe
    do i=is_pe,ie_pe
     k=max(kbot(i,j+1),kbot(i,j))
     if (k>0) diss(i,j,k) = maskV(i,j,k)*r_bot*v(i,j,k,tau)**2
    enddo
   enddo
   call calc_diss(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,diss,K_diss_bot,'V')
  endif
 endif

 if (.not.enable_hydrostatic) then
  !if (my_pe==0) print'(/a/)','ERROR: bottom friction for vertical velocity not implemented'
  !call halt_stop(' in bottom_friction')
 endif
end subroutine linear_bottom_friction





subroutine quadratic_bottom_friction
!=======================================================================
! quadratic bottom friction   
! dissipation is calculated and added to K_diss_bot
!=======================================================================
 use main_module   
 implicit none
 integer :: i,j,k
 real*8 :: diss(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz),fxa
 real*8 :: aloc(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx)


 ! we might want to account for EKE in the drag, also a tidal residual
 aloc=0.0
 do j=js_pe,je_pe
   do i=is_pe-1,ie_pe
    k=max(kbot(i,j),kbot(i+1,j))
    if (k>0) then
      fxa =       maskV(i  ,j,k)*v(i  ,j,k,tau)**2 + maskV(i  ,j-1,k)*v(i  ,j-1,k,tau)**2
      fxa = fxa + maskV(i+1,j,k)*v(i+1,j,k,tau)**2 + maskV(i+1,j-1,k)*v(i+1,j-1,k,tau)**2
      fxa = sqrt(u(i,j,k,tau)**2+ 0.25*fxa ) 
      aloc(i,j) =  maskU(i,j,k)*r_quad_bot*u(i,j,k,tau)*fxa/dzt(k)
      du_mix(i,j,k) = du_mix(i,j,k) - aloc(i,j)
    endif
   enddo
 enddo

 if (enable_conserve_energy) then
   diss=0.0
   do j=js_pe,je_pe
    do i=is_pe-1,ie_pe
     k=max(kbot(i,j),kbot(i+1,j))
     if (k>0) diss(i,j,k)= aloc(i,j)*u(i,j,k,tau)
    enddo
   enddo
   call calc_diss(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,diss,K_diss_bot,'U')
 endif
 
 aloc=0.0
 do j=js_pe-1,je_pe
   do i=is_pe,ie_pe
    k=max(kbot(i,j+1),kbot(i,j))
    if (k>0) then
      fxa =       maskU(i,j  ,k)*u(i,j  ,k,tau)**2 + maskU(i-1,j  ,k)*u(i-1,j  ,k,tau)**2
      fxa = fxa + maskU(i,j+1,k)*u(i,j+1,k,tau)**2 + maskU(i-1,j+1,k)*u(i-1,j+1,k,tau)**2
      fxa = sqrt(v(i,j,k,tau)**2+ 0.25*fxa ) 
      aloc(i,j)= maskV(i,j,k)*r_quad_bot*v(i,j,k,tau)*fxa/dzt(k)
      dv_mix(i,j,k)=dv_mix(i,j,k) - aloc(i,j)
    endif
   enddo
 enddo

 if (enable_conserve_energy) then
   diss=0.0
   do j=js_pe-1,je_pe
    do i=is_pe,ie_pe
     k=max(kbot(i,j+1),kbot(i,j))
     if (k>0) diss(i,j,k) = aloc(i,j)*v(i,j,k,tau)
    enddo
   enddo
   call calc_diss(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,diss,K_diss_bot,'V')
 endif

 if (.not.enable_hydrostatic) then
  !if (my_pe==0) print'(/a/)','ERROR: bottom friction for vertical velocity not implemented'
  !call halt_stop(' in quadratic_bottom_friction')
 endif
end subroutine quadratic_bottom_friction




subroutine momentum_sources_with_Euler
!=======================================================================
! other momentum sources with Euler forward time stepping
! dissipation is calculated and added to K_diss_bot
!=======================================================================
 use main_module   
 implicit none
 integer :: k
 real*8 :: diss(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz)

 do k=1,nz
  du_mix(:,:,k)=du_mix(:,:,k) + maskU(:,:,k)*u_source(:,:,k)
 enddo
 if (enable_conserve_energy) then
  do k=1,nz
   diss(:,:,k) = -maskU(:,:,k)*u(:,:,k,tau)*u_source(:,:,k)
  enddo
  call calc_diss(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,diss,K_diss_bot,'U')
 endif

 do k=1,nz
   dv_mix(:,:,k)=dv_mix(:,:,k) + maskV(:,:,k)*v_source(:,:,k)
 enddo
 if (enable_conserve_energy) then
  do k=1,nz
   diss(:,:,k) = -maskV(:,:,k)*v(:,:,k,tau)*v_source(:,:,k)
  enddo
  call calc_diss(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,diss,K_diss_bot,'V')
 endif
end subroutine momentum_sources_with_Euler



subroutine momentum_sources_with_AB
!=======================================================================
! other momentum sources with Adam-Bashforth time stepping
! dissipation is calculated and added to K_diss_bot
!=======================================================================
 use main_module   
 implicit none
 integer :: k
 real*8 :: diss(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz)

 do k=1,nz
  du(:,:,k,tau)=du(:,:,k,tau) + maskU(:,:,k)*u_source(:,:,k)
 enddo
 if (enable_conserve_energy) then
  do k=1,nz
   diss(:,:,k) = -maskU(:,:,k)*u(:,:,k,tau)*u_source(:,:,k)
  enddo
  call calc_diss(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,diss,K_diss_bot,'U')
 endif

 do k=1,nz
   dv(:,:,k,tau)=dv(:,:,k,tau) + maskV(:,:,k)*v_source(:,:,k)
 enddo
 if (enable_conserve_energy) then
  do k=1,nz
   diss(:,:,k) = -maskV(:,:,k)*v(:,:,k,tau)*v_source(:,:,k)
  enddo
  call calc_diss(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,diss,K_diss_bot,'V')
 endif
end subroutine momentum_sources_with_AB





