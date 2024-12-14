


subroutine momentum
!=======================================================================
! solve for momentum for taup1
!=======================================================================
 use main_module   
 use isoneutral_module   
 use idemix_module   
 use rossmix_module  
 use rossmix2_module  
 use timing_module   
 implicit none
 integer :: i,j,k

 !---------------------------------------------------------------------------------
 !  time tendency due to Coriolis force
 !---------------------------------------------------------------------------------
 do j=js_pe,je_pe
  do i=is_pe,ie_pe
   du_cor(i,j,:)= maskU(i,j,:)*(  coriolis_t(i  ,j)*(v(i  ,j,:,tau)+v(i  ,j-1,:,tau))*dxt(i  )/dxu(i)  &
                                 +coriolis_t(i+1,j)*(v(i+1,j,:,tau)+v(i+1,j-1,:,tau))*dxt(i+1)/dxu(i) )*0.25
   dv_cor(i,j,:)=-maskV(i,j,:)*(coriolis_t(i,j  )*(u(i-1,j  ,:,tau)+u(i,j  ,:,tau))*dyt(j  )*cost(j  )/( dyu(j)*cosu(j) )  &
                               +coriolis_t(i,j+1)*(u(i-1,j+1,:,tau)+u(i,j+1,:,tau))*dyt(j+1)*cost(j+1)/( dyu(j)*cosu(j) ) )*0.25
  enddo
 enddo
 !---------------------------------------------------------------------------------
 !  time tendency due to metric terms
 !---------------------------------------------------------------------------------
 if (coord_degree) then
  do j=js_pe,je_pe
   do i=is_pe,ie_pe
    du_cor(i,j,:) = du_cor(i,j,:) + maskU(i,j,:)*0.125*tantr(j)*( &
                      (u(i  ,j,:,tau)+u(i-1,j,:,tau))*(v(i  ,j,:,tau)+v(i  ,j-1,:,tau))*dxt(i  )/dxu(i) &
                    + (u(i+1,j,:,tau)+u(i  ,j,:,tau))*(v(i+1,j,:,tau)+v(i+1,j-1,:,tau))*dxt(i+1)/dxu(i) )
    dv_cor(i,j,:) = dv_cor(i,j,:) - maskV(i,j,:)*0.125*( &
                           tantr(j  )*(u(i,j  ,:,tau)+u(i-1,j  ,:,tau))**2*dyt(j  )*cost(j  )/( dyu(j)*cosu(j) ) &
                         + tantr(j+1)*(u(i,j+1,:,tau)+u(i-1,j+1,:,tau))**2*dyt(j+1)*cost(j+1)/( dyu(j)*cosu(j) ) )
   enddo
  enddo
 endif
 !---------------------------------------------------------------------------------
 ! non hydrostatic Coriolis terms, metric terms are neglected
 !---------------------------------------------------------------------------------
 if (.not. enable_hydrostatic) then
  do k=2,nz
   do i=is_pe,ie_pe
    du_cor(i,:,k) = du_cor(i,:,k) - maskU(i,:,k)*0.25*(coriolis_h(i  ,:)*area_t(i  ,:)*(w(i  ,:,k,tau)+w(i  ,:,k-1,tau)) &
                                                      +coriolis_h(i+1,:)*area_t(i+1,:)*(w(i+1,:,k,tau)+w(i+1,:,k-1,tau)) ) &
                                                        /area_u(i,:)
   enddo
  enddo
  k=1;
  do i=is_pe,ie_pe
   du_cor(i,:,k) = du_cor(i,:,k) - maskU(i,:,k)*0.25*(coriolis_h(i  ,:)*area_t(i  ,:)*(w(i  ,:,k,tau)) &
                                                     +coriolis_h(i+1,:)*area_t(i+1,:)*(w(i+1,:,k,tau)) )/area_u(i,:)
  enddo
  do k=1,nz-1
   do i=is_pe,ie_pe
    dw_cor(i,:,k) = maskW(i,:,k)*0.25*(coriolis_h(i,:)*dzt(k  )*(u(i,:,k  ,tau)+u(i-1,:,k  ,tau)) &
                                      +coriolis_h(i,:)*dzt(k+1)*(u(i,:,k+1,tau)+u(i-1,:,k+1,tau)) )/dzw(k)
   enddo
  enddo
 endif
 !---------------------------------------------------------------------------------
 ! transfer to time tendencies
 !---------------------------------------------------------------------------------
 do j=js_pe,je_pe
  do i=is_pe,ie_pe
   du(i,j,:,tau)=du_cor(i,j,:)
   dv(i,j,:,tau)=dv_cor(i,j,:)
  enddo
 enddo

 if (.not. enable_hydrostatic) then
  do j=js_pe,je_pe
   do i=is_pe,ie_pe
     dw(i,j,:,tau)=dw_cor(i,j,:)
   enddo
  enddo
 endif

 !---------------------------------------------------------------------------------
 ! wind stress forcing
 !---------------------------------------------------------------------------------
 do j=js_pe,je_pe
  do i=is_pe,ie_pe
   du(i,j,nz,tau)=du(i,j,nz,tau)+maskU(i,j,nz)*surface_taux(i,j)/dzt(nz)
   dv(i,j,nz,tau)=dv(i,j,nz,tau)+maskV(i,j,nz)*surface_tauy(i,j)/dzt(nz)
  enddo
 enddo

 !---------------------------------------------------------------------------------
 ! advection
 !---------------------------------------------------------------------------------
 call momentum_advection
 du(:,:,:,tau) = du(:,:,:,tau) + du_adv
 dv(:,:,:,tau) = dv(:,:,:,tau) + dv_adv
 if (.not. enable_hydrostatic) dw(:,:,:,tau) = dw(:,:,:,tau) + dw_adv

 call tic('fric')
 !---------------------------------------------------------------------------------
 ! vertical friction
 !---------------------------------------------------------------------------------
 K_diss_v = 0d0
 if (enable_implicit_vert_friction) call implicit_vert_friction
 if (enable_explicit_vert_friction) call explicit_vert_friction

 !---------------------------------------------------------------------------------
 ! TEM formalism for eddy-driven velocity
 !---------------------------------------------------------------------------------
 if (enable_TEM_friction) call isoneutral_friction

 !---------------------------------------------------------------------------------
 ! vertical friction by internal waves
 !---------------------------------------------------------------------------------
 if (.not. enable_implicit_vert_friction.and..not.enable_explicit_vert_friction) then
   du_mix=0; dv_mix=0
 endif
 if (enable_idemix .and. enable_idemix3)  call idemix_friction(tau)
 if (enable_idemix .and. enable_leewaves) call idemix_friction_lee
 
 if (enable_rossmix .and. enable_rossmix_mean_flow_interaction ) call rossmix_friction
 if (enable_rossmix2  ) call rossmix2_friction
 !---------------------------------------------------------------------------------
 !horizontal friction
 !---------------------------------------------------------------------------------
 if (enable_conserve_energy) K_diss_h = 0d0
 if (enable_hor_friction)        call harmonic_friction
 if (enable_biharmonic_friction) call biharmonic_friction
 call biharmonic_thickness_friction
 
 !---------------------------------------------------------------------------------
 ! Rayleigh and bottom friction
 !---------------------------------------------------------------------------------
 if (enable_conserve_energy) K_diss_bot = 0d0
 if (enable_ray_friction)              call rayleigh_friction
 if (enable_bottom_friction)           call linear_bottom_friction
 if (enable_quadratic_bottom_friction) call quadratic_bottom_friction

 !---------------------------------------------------------------------------------
 ! add user defined forcing
 !---------------------------------------------------------------------------------
 if (enable_momentum_sources)         call momentum_sources_with_Euler
 if (enable_momentum_sources_with_AB) call momentum_sources_with_AB
 call toc('fric')

 !---------------------------------------------------------------------------------
 ! account for open boundaries
 !---------------------------------------------------------------------------------
 call set_obc_momentum
 
 !---------------------------------------------------------------------------------
 ! external mode
 !---------------------------------------------------------------------------------
 call tic('press')
 if (enable_streamfunction) then
  call solve_streamfunction
 else
  call solve_pressure
  if (itt==0) then 
       psi(:,:,tau)=psi(:,:,taup1)
       psi(:,:,taum1)=psi(:,:,taup1)
  endif
 endif
 if (.not. enable_hydrostatic) call solve_non_hydrostatic
 call toc('press')

end subroutine momentum



