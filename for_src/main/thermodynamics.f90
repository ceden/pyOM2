




subroutine thermodynamics
!=======================================================================
! integrate temperature and salinity and diagnose sources of dynamic enthalpy
!=======================================================================
 use main_module   
 use isoneutral_module   
 use tke_module   
 use rossmix_module
 use rossmix2_module
 use timing_module   
 implicit none
 integer :: i,j,k,ks
 real*8 :: a_tri(nz),b_tri(nz),c_tri(nz),d_tri(nz),delta(nz)
 real*8 :: aloc(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz),fxa,fxb
 real*8 :: get_drhodT,get_drhodS

 call advect_tracer(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,temp(:,:,:,tau),dtemp(:,:,:,tau))
 call advect_tracer(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,salt(:,:,:,tau),dsalt(:,:,:,tau))
 
 !call advect_temperature
 !call advect_salinity

 if (enable_conserve_energy) call  calc_dissipation_by_advection

 if (enable_conserve_energy .and. enable_tke .and. enable_take_p_diss_adv_from_tke) then
 !---------------------------------------------------------------------------------
 ! distribute P_diss_adv over domain, prevent draining of TKE
 !---------------------------------------------------------------------------------
  fxa =0; fxb = 0
  do k=1,nz-1
   do j=js_pe,je_pe
    do i=is_pe,ie_pe
     fxa = fxa + area_t(i,j)*P_diss_adv(i,j,k)*dzw(k)*maskW(i,j,k)
     if (tke(i,j,k,tau) > 0.0 ) fxb = fxb + area_t(i,j)*dzw(k)*maskW(i,j,k) 
    enddo
   enddo
  enddo
  k=nz
  do j=js_pe,je_pe
   do i=is_pe,ie_pe
    fxa = fxa + 0.5*area_t(i,j)*P_diss_adv(i,j,k)*dzw(k)*maskW(i,j,k)
    fxb = fxb + 0.5*area_t(i,j)*dzw(k)*maskW(i,j,k)
   enddo
  enddo
  call global_sum(fxa); call global_sum(fxb)
  P_diss_adv = 0.0
  do k=1,nz 
   do j=js_pe,je_pe
    do i=is_pe,ie_pe
     if (tke(i,j,k,tau) > 0.0 .or. k==nz) P_diss_adv(i,j,k) = fxa/fxb
    enddo
   enddo
  enddo
 endif

 if (enable_AB_time_stepping) then
   !---------------------------------------------------------------------------------
   ! Adam Bashforth time stepping for advection
   !---------------------------------------------------------------------------------
   temp(:,:,:,taup1)=temp(:,:,:,tau)+dt_tracer*( (1.5+AB_eps)*dtemp(:,:,:,tau) - ( 0.5+AB_eps)*dtemp(:,:,:,taum1))*maskT
   salt(:,:,:,taup1)=salt(:,:,:,tau)+dt_tracer*( (1.5+AB_eps)*dsalt(:,:,:,tau) - ( 0.5+AB_eps)*dsalt(:,:,:,taum1))*maskT
 else
   temp(:,:,:,taup1)=temp(:,:,:,tau)+dt_tracer*dtemp(:,:,:,tau)*maskT
   salt(:,:,:,taup1)=salt(:,:,:,tau)+dt_tracer*dsalt(:,:,:,tau)*maskT 
 endif

 !---------------------------------------------------------------------------------
 ! horizontal diffusion
 !---------------------------------------------------------------------------------
 call tic('iso')
 if  (enable_hor_diffusion.or.enable_biharmonic_mixing) then
    dtemp_hmix = 0d0; dsalt_hmix = 0d0
    if (enable_conserve_energy) P_diss_hmix = 0d0
 endif
 if (enable_hor_diffusion)     call tempsalt_diffusion
 if (enable_biharmonic_mixing) call tempsalt_biharmonic

 !---------------------------------------------------------------------------------
 ! sources like restoring zones, etc
 !---------------------------------------------------------------------------------
 if (enable_tempsalt_sources) call tempsalt_sources

 !---------------------------------------------------------------------------------
 ! isopycnal diffusion
 !---------------------------------------------------------------------------------
 if (enable_neutral_diffusion) then
  P_diss_iso = 0.0; dtemp_iso = 0.0; dsalt_iso = 0.0
  call isoneutral_diffusion_pre
  call isoneutral_diffusion(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,temp,.true.)
  call isoneutral_diffusion(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,salt,.false.)
  if (enable_skew_diffusion) then
   P_diss_skew = 0.0;
   call isoneutral_skew_diffusion(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,temp,.true.)
   call isoneutral_skew_diffusion(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,salt,.false.)
  endif
 endif

 !---------------------------------------------------------------------------------
 ! mixed layer eddy parameterization
 !---------------------------------------------------------------------------------
 if (enable_ml_para) call mixed_layer_para

 !---------------------------------------------------------------------------------
 !---------------------------------------------------------------------------------
 !if (enable_rossmix .and. enable_rossmix_bolus_form) call rossmix_eddy_advect
 if (enable_rossmix) call rossmix_eddy_advect
 if (enable_rossmix2) call rossmix2_eddy_advect
 call biharmonic_thickness_mixing
 
 
 call toc('iso')

 call tic('vmix')
 !---------------------------------------------------------------------------------
 ! vertical mixing of temperature and salinity
 !---------------------------------------------------------------------------------
 dtemp_vmix = temp(:,:,:,taup1) ; dsalt_vmix = salt(:,:,:,taup1) 
 a_tri=0.0;b_tri=0.0; c_tri=0.0; d_tri=0.0; delta=0.0
 do j=js_pe,je_pe
   do i=is_pe,ie_pe
    ks=kbot(i,j)
    if (ks>0) then
     do k=ks,nz-1
      delta(k) = dt_tracer/dzw(k)*kappaH(i,j,k)
     enddo
     delta(nz)=0.0
     do k=ks+1,nz
       a_tri(k) = - delta(k-1)/dzt(k)
     enddo
     a_tri(ks)=0.0
     do k=ks+1,nz-1
      b_tri(k) = 1+ delta(k)/dzt(k) + delta(k-1)/dzt(k) 
     enddo
     b_tri(nz) = 1+ delta(nz-1)/dzt(nz) 
     b_tri(ks) = 1+ delta(ks)/dzt(ks)   
     do k=ks,nz-1
      c_tri(k) = - delta(k)/dzt(k)
     enddo
     c_tri(nz)=0.0
     d_tri(ks:nz)=temp(i,j,ks:nz,taup1) 
     d_tri(nz) = d_tri(nz) + dt_tracer*forc_temp_surface(i,j)/dzt(nz)
     call solve_tridiag(a_tri(ks:nz),b_tri(ks:nz),c_tri(ks:nz),d_tri(ks:nz),temp(i,j,ks:nz,taup1),nz-ks+1)
     d_tri(ks:nz)=salt(i,j,ks:nz,taup1) 
     d_tri(nz) = d_tri(nz) + dt_tracer*forc_salt_surface(i,j)/dzt(nz)
     call solve_tridiag(a_tri(ks:nz),b_tri(ks:nz),c_tri(ks:nz),d_tri(ks:nz),salt(i,j,ks:nz,taup1),nz-ks+1)
    endif
   enddo
 enddo
 dtemp_vmix = (temp(:,:,:,taup1)-dtemp_vmix)/dt_tracer; dsalt_vmix = (salt(:,:,:,taup1)-dsalt_vmix)/dt_tracer
 call toc('vmix')

 !---------------------------------------------------------------------------------
 ! open boundaries
 !---------------------------------------------------------------------------------
 call set_obc_tempsalt
 
 !---------------------------------------------------------------------------------
 ! boundary exchange
 !---------------------------------------------------------------------------------
 call border_exchg_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,temp(:,:,:,taup1)) 
 call setcyclic_xyz   (is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,temp(:,:,:,taup1))
 call border_exchg_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,salt(:,:,:,taup1)) 
 call setcyclic_xyz   (is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,salt(:,:,:,taup1))
 call set_obc_boundary_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,temp(:,:,:,taup1)) 
 call set_obc_boundary_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,salt(:,:,:,taup1)) 

 call tic('eq_of_state')
 call calc_eq_of_state(taup1)
 call toc('eq_of_state')

 !---------------------------------------------------------------------------------
 ! surface density flux
 !---------------------------------------------------------------------------------
 do j=js_pe-onx,je_pe+onx
  do i=is_pe-onx,ie_pe+onx
    forc_rho_surface(i,j)=(get_drhodT(salt(i,j,nz,taup1),temp(i,j,nz,taup1),abs(zt(nz)))*forc_temp_surface(i,j) & 
                          +get_drhodS(salt(i,j,nz,taup1),temp(i,j,nz,taup1),abs(zt(nz)))*forc_salt_surface(i,j) )  &
                                    *maskT(i,j,nz)
  enddo
 enddo

 call tic('vmix')
 P_diss_v = 0.0 
 if (enable_conserve_energy) then
 !---------------------------------------------------------------------------------
 ! diagnose dissipation of dynamic enthalpy by vertical mixing
 !---------------------------------------------------------------------------------
  do k=1,nz-1
   do j=js_pe,je_pe
    do i=is_pe,ie_pe
     fxa = (-int_drhodT(i,j,k+1,taup1) +int_drhodT(i,j,k,taup1))/dzw(k)
     P_diss_v(i,j,k)=P_diss_v(i,j,k) &
                    -grav/rho_0*fxa*kappaH(i,j,k)*(temp(i,j,k+1,taup1)-temp(i,j,k,taup1))/dzw(k)*maskW(i,j,k)
     fxa = (-int_drhodS(i,j,k+1,taup1) +int_drhodS(i,j,k,taup1))/dzw(k)
     P_diss_v(i,j,k)= P_diss_v(i,j,k) &
                    -grav/rho_0*fxa*kappaH(i,j,k)*(salt(i,j,k+1,taup1)-salt(i,j,k,taup1))/dzw(k)*maskW(i,j,k)
    enddo
   enddo
  end do
  k=nz
  do j=js_pe,je_pe
   do i=is_pe,ie_pe
    fxa =  2*int_drhodT(i,j,k,taup1)/dzw(k)
    P_diss_v(i,j,k)=P_diss_v(i,j,k)  -grav/rho_0*fxa*forc_temp_surface(i,j)*maskW(i,j,k)
    fxa =  2*int_drhodS(i,j,k,taup1)/dzw(k)
    P_diss_v(i,j,k)= P_diss_v(i,j,k) -grav/rho_0*fxa*forc_salt_surface(i,j)*maskW(i,j,k)
   enddo
  enddo
 endif

 if (enable_conserve_energy) then
 !---------------------------------------------------------------------------------
 ! determine effect due to nonlinear equation of state
 !---------------------------------------------------------------------------------
  aloc(:,:,1:nz-1)=kappaH(:,:,1:nz-1)*Nsqr(:,:,1:nz-1,taup1)  
  P_diss_nonlin(:,:,1:nz-1) = P_diss_v(:,:,1:nz-1)-aloc(:,:,1:nz-1)
  P_diss_v(:,:,1:nz-1)      = aloc(:,:,1:nz-1)
 else
 !---------------------------------------------------------------------------------
 ! diagnose N^2 kappaH, i.e. exchange of pot. energy with TKE
 !---------------------------------------------------------------------------------
  P_diss_v(:,:,1:nz-1) = kappaH(:,:,1:nz-1)*Nsqr(:,:,1:nz-1,taup1)  
  P_diss_v(:,:,nz)=-forc_rho_surface(:,:)*maskT(:,:,nz)*grav/rho_0
 endif
 call toc('vmix')

end subroutine thermodynamics






subroutine calc_dissipation_by_advection
 !---------------------------------------------------------------------------------
 ! diagnose sources of dynamic enthalpy by advection of T and S
 !---------------------------------------------------------------------------------
 use main_module   
 implicit none
 integer :: i,j,k,ks
 real*8 :: aloc(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz),fxa

 !---------------------------------------------------------------------------------
 ! advection of dynamic enthalpy
 !---------------------------------------------------------------------------------
 call advect_tracer(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,Hd(:,:,:,tau),dHd(:,:,:,tau))

 !---------------------------------------------------------------------------------
 ! changes in dyn. Enthalpy due to advection
 !---------------------------------------------------------------------------------
  do k=1,nz
   do j=js_pe,je_pe
    do i=is_pe,ie_pe  
     fxa =  grav/rho_0*( -int_drhodT(i,j,k,tau)*dtemp(i,j,k,tau) -int_drhodS(i,j,k,tau)*dsalt(i,j,k,tau)  ) 
     aloc(i,j,k) =fxa - dHd(i,j,k,tau) 
    enddo
   enddo
  enddo
 !---------------------------------------------------------------------------------
 ! contribution by vertical advection is - g rho w /rho0, substract this also
 !---------------------------------------------------------------------------------
  do k=1,nz-1
   aloc(:,:,k) = aloc(:,:,k) - 0.25*grav/rho_0*w(:,:,k  ,tau)*(rho(:,:,k,tau)+rho(:,:,k+1,tau))*dzw(k)/dzt(k)  
  enddo
  do k=2,nz
   aloc(:,:,k) = aloc(:,:,k) - 0.25*grav/rho_0*w(:,:,k-1,tau)*(rho(:,:,k,tau)+rho(:,:,k-1,tau))*dzw(k-1)/dzt(k)
  enddo

 !---------------------------------------------------------------------------------
 ! dissipation by advection interpolated on W-grid
 !---------------------------------------------------------------------------------
  do j=js_pe,je_pe
   do i=is_pe,ie_pe
     ks=kbot(i,j)
     if (ks>0) then
      k=ks; P_diss_adv(i,j,k) = 0.5*(aloc(i,j,k)+aloc(i,j,k+1)) + 0.5*aloc(i,j,k)*dzw(max(1,k-1))/dzw(k)
      do k=ks+1,nz-1
       P_diss_adv(i,j,k) =  0.5*(aloc(i,j,k) +aloc(i,j,k+1))
      enddo
      k=nz; P_diss_adv(i,j,k) = aloc(i,j,k)
     endif
   enddo
  enddo


end subroutine calc_dissipation_by_advection



subroutine calc_eq_of_state(n)
!=======================================================================
! calculate density, stability frequency, dynamic enthalpy and derivatives 
! for time level n from temperature and salinity
!=======================================================================
 use main_module   
 implicit none
 integer, intent(in) :: n
 integer :: i,j,k
 real*8 :: get_rho,get_int_drhodT,get_int_drhodS, get_dyn_enthalpy
 real*8 :: fxa

 !---------------------------------------------------------------------------------
 ! calculate new density
 !---------------------------------------------------------------------------------
 do k=1,nz
  do j=js_pe-onx,je_pe+onx
   do i=is_pe-onx,ie_pe+onx
    rho(i,j,k,n) = get_rho(salt(i,j,k,n),temp(i,j,k,n),abs(zt(k)))*maskT(i,j,k)
   enddo
  enddo
 enddo

 if (enable_conserve_energy) then
 !---------------------------------------------------------------------------------
 ! calculate new dynamic enthalpy and derivatives
 !---------------------------------------------------------------------------------
  do k=1,nz
   do j=js_pe-onx,je_pe+onx
    do i=is_pe-onx,ie_pe+onx
     Hd(i,j,k,n) = get_dyn_enthalpy(salt(i,j,k,n),temp(i,j,k,n),abs(zt(k)))*maskT(i,j,k)
     int_drhodT(i,j,k,n) = get_int_drhodT(salt(i,j,k,n),temp(i,j,k,n),abs(zt(k)))
     int_drhodS(i,j,k,n) = get_int_drhodS(salt(i,j,k,n),temp(i,j,k,n),abs(zt(k)))
    enddo
   enddo
  enddo
 endif

 !---------------------------------------------------------------------------------
 ! new stability frequency
 !---------------------------------------------------------------------------------
 do k=1,nz-1
  do j=js_pe-onx,je_pe+onx
   do i=is_pe-onx,ie_pe+onx
     fxa =  -grav/rho_0/dzw(k)*maskW(i,j,k)
     Nsqr(i,j,k,n) =fxa*(get_rho(salt(i,j,k+1,n),temp(i,j,k+1,n),abs(zt(k)))-rho(i,j,k,n))
   enddo
  enddo
 enddo
 Nsqr(:,:,nz,n)=Nsqr(:,:,nz-1,n)

end subroutine calc_eq_of_state

