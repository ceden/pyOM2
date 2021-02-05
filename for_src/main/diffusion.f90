



subroutine tempsalt_biharmonic
!---------------------------------------------------------------------------------
! biharmonic mixing of temp and salinity, 
! dissipation of dyn. Enthalpy is stored
!---------------------------------------------------------------------------------
 use main_module   
 implicit none
 integer :: i,j,k,ks,is,ie,js,je
 real*8 :: aloc(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz),diff
 real*8 :: del2(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz),fxa

 is = is_pe-onx; ie = ie_pe+onx; js = js_pe-onx; je = je_pe+onx
 diff = sqrt(abs(K_hbi))

 !---------------------------------------------------------------------------------
 ! Temperature
 !---------------------------------------------------------------------------------
 do j=js,je
   fxa = cost(j)**biha_mixing_cosPower
   do i=is,ie-1
    flux_east(i,j,:)=-diff*fxa*(temp(i+1,j,:,tau)-temp(i,j,:,tau))/(cost(j)*dxu(i))*maskU(i,j,:)
   enddo
 enddo
 do j=js,je-1
    fxa = cosu(j)**biha_mixing_cosPower
    flux_north(:,j,:)=-diff*fxa*(temp(:,j+1,:,tau)-temp(:,j,:,tau))/dyu(j)*maskV(:,j,:)*cosu(j)
 enddo 
 flux_east(ie,:,:)=0.; flux_north(:,je,:)=0.

 do j=js+1,je
   do i=is+1,ie
    del2(i,j,:)= maskT(i,j,:)* (flux_east(i,j,:) - flux_east(i-1,j,:))/(cost(j)*dxt(i)) &
                              +(flux_north(i,j,:) - flux_north(i,j-1,:))/(cost(j)*dyt(j)) 
  enddo
 enddo
 call border_exchg_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,del2) 
 call setcyclic_xyz   (is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,del2)

 do j=js,je
   fxa = cost(j)**biha_mixing_cosPower
   do i=is,ie-1
    flux_east(i,j,:)=diff*fxa*(del2(i+1,j,:)-del2(i,j,:))/(cost(j)*dxu(i))*maskU(i,j,:)
   enddo
 enddo
 do j=js,je-1
    fxa = cosu(j)**biha_mixing_cosPower
    flux_north(:,j,:)=diff*fxa*(del2(:,j+1,:)-del2(:,j,:))/dyu(j)*maskV(:,j,:)*cosu(j)
 enddo 
 flux_east(ie,:,:)=0.; flux_north(:,je,:)=0.

 ! update tendency
 do k=1,nz
  do j=js_pe,je_pe
   do i=is_pe,ie_pe
      fxa= maskT(i,j,k)*( (flux_east(i,j,k) - flux_east(i-1,j,k))/(cost(j)*dxt(i)) &
                         +(flux_north(i,j,k) - flux_north(i,j-1,k))/(cost(j)*dyt(j))  )
      dtemp_hmix(i,j,k)=  dtemp_hmix(i,j,k) + fxa
      temp(i,j,k,taup1)=temp(i,j,k,taup1)+dt_tracer*fxa
   enddo
  enddo
 enddo

if (enable_conserve_energy) then
  ! diagnose dissipation of dynamic enthalpy by hor. mixing of temperature
  do k=1,nz
   do j=js_pe,je_pe
    do i=is_pe,ie_pe
     fxa = int_drhodT(i,j,k,tau)
     aloc(i,j,k) =+0.5*grav/rho_0*( (int_drhodT(i+1,j,k,tau)-fxa)*flux_east(i  ,j,k) &
                                   +(fxa-int_drhodT(i-1,j,k,tau))*flux_east(i-1,j,k) ) /(dxt(i)*cost(j))  &
                  +0.5*grav/rho_0*( (int_drhodT(i,j+1,k,tau)-fxa)*flux_north(i,j  ,k) &
                                   +(fxa-int_drhodT(i,j-1,k,tau))*flux_north(i,j-1,k) ) /(dyt(j)*cost(j)) 
    enddo
   enddo
  end do
  ! dissipation interpolated on W-grid
  do j=js_pe,je_pe
   do i=is_pe,ie_pe
     ks=kbot(i,j)
     if (ks>0) then
        k=ks; P_diss_hmix(i,j,k) =  P_diss_hmix(i,j,k) + &
                   0.5*(aloc(i,j,k)+aloc(i,j,k+1)) + 0.5*aloc(i,j,k)*dzw(max(1,k-1))/dzw(k)
      do k=ks+1,nz-1
       P_diss_hmix(i,j,k) = P_diss_hmix(i,j,k) + 0.5*(aloc(i,j,k) +aloc(i,j,k+1))
      enddo
      k=nz; P_diss_hmix(i,j,k) =  P_diss_hmix(i,j,k) + aloc(i,j,k)
     endif
   enddo
  enddo
endif

 !---------------------------------------------------------------------------------
 ! Salinity
 !---------------------------------------------------------------------------------
 do j=js,je
   fxa = cost(j)**biha_mixing_cosPower
   do i=is,ie-1
    flux_east(i,j,:)=-diff*fxa*(salt(i+1,j,:,tau)-salt(i,j,:,tau))/(cost(j)*dxu(i))*maskU(i,j,:)
   enddo
 enddo
 do j=js,je-1
    fxa = cosu(j)**biha_mixing_cosPower
    flux_north(:,j,:)=-diff*fxa*(salt(:,j+1,:,tau)-salt(:,j,:,tau))/dyu(j)*maskV(:,j,:)*cosu(j)
 enddo 
 flux_east(ie,:,:)=0.; flux_north(:,je,:)=0.

 do j=js+1,je
   do i=is+1,ie
    del2(i,j,:)= maskT(i,j,:)* (flux_east(i,j,:) - flux_east(i-1,j,:))/(cost(j)*dxt(i)) &
                              +(flux_north(i,j,:) - flux_north(i,j-1,:))/(cost(j)*dyt(j)) 
  enddo
 enddo
 call border_exchg_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,del2) 
 call setcyclic_xyz   (is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,del2)

 do j=js,je
   fxa = cost(j)**biha_mixing_cosPower
   do i=is,ie-1
    flux_east(i,j,:)=diff*fxa*(del2(i+1,j,:)-del2(i,j,:))/(cost(j)*dxu(i))*maskU(i,j,:)
   enddo
 enddo
 do j=js,je-1
    fxa = cosu(j)**biha_mixing_cosPower
    flux_north(:,j,:)=diff*fxa*(del2(:,j+1,:)-del2(:,j,:))/dyu(j)*maskV(:,j,:)*cosu(j)
 enddo 
 flux_east(ie,:,:)=0.; flux_north(:,je,:)=0.

 ! update tendency 
 do k=1,nz
  do j=js_pe,je_pe
   do i=is_pe,ie_pe
      fxa= maskT(i,j,k)*( (flux_east(i,j,k) - flux_east(i-1,j,k))/(cost(j)*dxt(i)) &
                         +(flux_north(i,j,k) - flux_north(i,j-1,k))/(cost(j)*dyt(j))  )
      dsalt_hmix(i,j,k)=  dsalt_hmix(i,j,k) + fxa
      salt(i,j,k,taup1)=salt(i,j,k,taup1)+dt_tracer*fxa
   enddo
  enddo
 enddo


if (enable_conserve_energy) then
  ! diagnose dissipation of dynamic enthalpy by hor. mixing of salinity
  do k=1,nz
   do j=js_pe-onx+1,je_pe+onx-1
    do i=is_pe-onx+1,ie_pe+onx-1
     fxa = int_drhodS(i,j,k,tau)
     aloc(i,j,k) =+0.5*grav/rho_0*( (int_drhodS(i+1,j,k,tau)-fxa)*flux_east(i  ,j,k) &
                                   +(fxa-int_drhodS(i-1,j,k,tau))*flux_east(i-1,j,k) ) /(dxt(i)*cost(j))  &
                  +0.5*grav/rho_0*( (int_drhodS(i,j+1,k,tau)-fxa)*flux_north(i,j  ,k) &
                                   +(fxa-int_drhodS(i,j-1,k,tau))*flux_north(i,j-1,k) ) /(dyt(j)*cost(j)) 
    enddo
   enddo
  end do
  ! dissipation interpolated on W-grid
  do j=js_pe-onx,je_pe+onx
   do i=is_pe-onx,ie_pe+onx
     ks=kbot(i,j)
     if (ks>0) then
      k=ks; P_diss_hmix(i,j,k) = P_diss_hmix(i,j,k)+ &
                          0.5*(aloc(i,j,k)+aloc(i,j,k+1)) + 0.5*aloc(i,j,k)*dzw(max(1,k-1))/dzw(k)
      do k=ks+1,nz-1
       P_diss_hmix(i,j,k) =  P_diss_hmix(i,j,k)+ 0.5*(aloc(i,j,k) +aloc(i,j,k+1))
      enddo
      k=nz; P_diss_hmix(i,j,k) = P_diss_hmix(i,j,k)+ aloc(i,j,k)
     endif
   enddo
  enddo
endif
end subroutine tempsalt_biharmonic








subroutine tempsalt_diffusion
!---------------------------------------------------------------------------------
! Diffusion of temp and salinity, 
! dissipation of dyn. Enthalpy is stored
!---------------------------------------------------------------------------------
 use main_module   
 implicit none
 integer :: i,j,k,ks
 real*8 :: fxa
 real*8 :: aloc(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz)

  ! horizontal diffusion of temperature
  do j=js_pe-onx,je_pe+onx
   do i=is_pe-onx,ie_pe+onx-1
    flux_east(i,j,:)=K_h*(temp(i+1,j,:,tau)-temp(i,j,:,tau))/(cost(j)*dxu(i))*maskU(i,j,:)
   enddo
  enddo
  flux_east(ie_pe+onx,:,:)=0.
  do j=js_pe-onx,je_pe+onx-1
    flux_north(:,j,:)=K_h*(temp(:,j+1,:,tau)-temp(:,j,:,tau))/dyu(j)*maskV(:,j,:)*cosu(j)
  enddo
  flux_north(:,je_pe+onx,:)=0.

  if (enable_hor_friction_cos_scaling) then
   do j=js_pe-onx,je_pe+onx
    flux_east(:,j,:)  = flux_east(:,j,:)*cost(j)**hor_friction_cosPower
    flux_north(:,j,:) = flux_north(:,j,:)*cosu(j)**hor_friction_cosPower
   enddo
  endif

  do k=1,nz
   do j=js_pe,je_pe
    do i=is_pe,ie_pe
      fxa = maskT(i,j,k)*((flux_east(i,j,k) - flux_east(i-1,j,k))/(cost(j)*dxt(i)) &
                         +(flux_north(i,j,k) -flux_north(i,j-1,k))/(cost(j)*dyt(j)) )
      dtemp_hmix(i,j,k)= dtemp_hmix(i,j,k)+ fxa
      temp(i,j,k,taup1)=temp(i,j,k,taup1)+dt_tracer*fxa
    enddo
   enddo
  enddo


if (enable_conserve_energy) then
  ! diagnose dissipation of dynamic enthalpy by hor. mixing of temperature
  do k=1,nz
   do j=js_pe-onx+1,je_pe+onx-1
    do i=is_pe-onx+1,ie_pe+onx-1
     fxa = int_drhodT(i,j,k,tau)
     aloc(i,j,k) =+0.5*grav/rho_0*( (int_drhodT(i+1,j,k,tau)-fxa)*flux_east(i  ,j,k) &
                                   +(fxa-int_drhodT(i-1,j,k,tau))*flux_east(i-1,j,k) ) /(dxt(i)*cost(j))  &
                  +0.5*grav/rho_0*( (int_drhodT(i,j+1,k,tau)-fxa)*flux_north(i,j  ,k) &
                                   +(fxa-int_drhodT(i,j-1,k,tau))*flux_north(i,j-1,k) ) /(dyt(j)*cost(j)) 
    enddo
   enddo
  end do
  ! dissipation interpolated on W-grid
  do j=js_pe-onx,je_pe+onx
   do i=is_pe-onx,ie_pe+onx
     ks=kbot(i,j)
     if (ks>0) then
        k=ks; P_diss_hmix(i,j,k) =  P_diss_hmix(i,j,k)  + &
              0.5*(aloc(i,j,k)+aloc(i,j,k+1)) + 0.5*aloc(i,j,k)*dzw(max(1,k-1))/dzw(k)
      do k=ks+1,nz-1
       P_diss_hmix(i,j,k) = P_diss_hmix(i,j,k) + 0.5*(aloc(i,j,k) +aloc(i,j,k+1))
      enddo
      k=nz; P_diss_hmix(i,j,k) =  P_diss_hmix(i,j,k) + aloc(i,j,k)
     endif
   enddo
  enddo
endif

  ! horizontal diffusion of salinity
  do j=js_pe-onx,je_pe+onx
   do i=is_pe-onx,ie_pe+onx-1
    flux_east(i,j,:)=K_h*(salt(i+1,j,:,tau)-salt(i,j,:,tau))/(cost(j)*dxu(i))*maskU(i,j,:)
   enddo
  enddo
  flux_east(ie_pe+onx,:,:)=0.
  do j=js_pe-onx,je_pe+onx-1
    flux_north(:,j,:)=K_h*(salt(:,j+1,:,tau)-salt(:,j,:,tau))/dyu(j)*maskV(:,j,:)*cosu(j)
  enddo
  flux_north(:,je_pe+onx,:)=0.

  if (enable_hor_friction_cos_scaling) then
   do j=js_pe-onx,je_pe+onx
    flux_east(:,j,:)  = flux_east(:,j,:)*cost(j)**hor_friction_cosPower
    flux_north(:,j,:) = flux_north(:,j,:)*cosu(j)**hor_friction_cosPower
   enddo
  endif

  do k=1,nz
   do j=js_pe,je_pe
    do i=is_pe,ie_pe
      fxa= maskT(i,j,k)*( (flux_east(i,j,k) - flux_east(i-1,j,k))/(cost(j)*dxt(i)) &
                         +(flux_north(i,j,k) -flux_north(i,j-1,k))/(cost(j)*dyt(j)) )
      dsalt_hmix(i,j,k)= dsalt_hmix(i,j,k) + fxa
      salt(i,j,k,taup1)=salt(i,j,k,taup1)+dt_tracer*fxa
    enddo
   enddo
  enddo


if (enable_conserve_energy) then
  ! diagnose dissipation of dynamic enthalpy by hor. mixing of salinity
  do k=1,nz
   do j=js_pe-onx+1,je_pe+onx-1
    do i=is_pe-onx+1,ie_pe+onx-1
     fxa = int_drhodS(i,j,k,tau)
     aloc(i,j,k) =+0.5*grav/rho_0*( (int_drhodS(i+1,j,k,tau)-fxa)*flux_east(i  ,j,k) &
                                   +(fxa-int_drhodS(i-1,j,k,tau))*flux_east(i-1,j,k) ) /(dxt(i)*cost(j))  &
                  +0.5*grav/rho_0*( (int_drhodS(i,j+1,k,tau)-fxa)*flux_north(i,j  ,k) &
                                   +(fxa-int_drhodS(i,j-1,k,tau))*flux_north(i,j-1,k) ) /(dyt(j)*cost(j)) 
    enddo
   enddo
  end do
  ! dissipation interpolated on W-grid
  do j=js_pe-onx,je_pe+onx
   do i=is_pe-onx,ie_pe+onx
     ks=kbot(i,j)
     if (ks>0) then
      k=ks; P_diss_hmix(i,j,k) = P_diss_hmix(i,j,k)+ &
                          0.5*(aloc(i,j,k)+aloc(i,j,k+1)) + 0.5*aloc(i,j,k)*dzw(max(1,k-1))/dzw(k)
      do k=ks+1,nz-1
       P_diss_hmix(i,j,k) =  P_diss_hmix(i,j,k)+ 0.5*(aloc(i,j,k) +aloc(i,j,k+1))
      enddo
      k=nz; P_diss_hmix(i,j,k) = P_diss_hmix(i,j,k)+ aloc(i,j,k)
     endif
   enddo
  enddo
endif
end subroutine tempsalt_diffusion






subroutine tracer_biharmonic(is_,ie_,js_,je_,nz_,tra)
!---------------------------------------------------------------------------------
! biharmonic mixing of tracer, taup1 is only updated but not set
!---------------------------------------------------------------------------------
 use main_module   
 implicit none
 integer :: is_,ie_,js_,je_,nz_
 real*8, dimension(is_:ie_,js_:je_,nz_,3) :: tra
 integer :: i,j,k,is,ie,js,je
 real*8 :: del2(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz),fxa,diff

 is = is_pe-onx; ie = ie_pe+onx; js = js_pe-onx; je = je_pe+onx
 diff = sqrt(abs(K_hbi))

 do j=js,je
   fxa = cost(j)**biha_mixing_cosPower
   do i=is,ie-1
    flux_east(i,j,:)=-diff*fxa*(tra(i+1,j,:,tau)-tra(i,j,:,tau))/(cost(j)*dxu(i))*maskU(i,j,:)
   enddo
 enddo
 do j=js,je-1
    fxa = cosu(j)**biha_mixing_cosPower
    flux_north(:,j,:)=-diff*fxa*(tra(:,j+1,:,tau)-tra(:,j,:,tau))/dyu(j)*maskV(:,j,:)*cosu(j)
 enddo 
 flux_east(ie,:,:)=0.; flux_north(:,je,:)=0.

 do j=js+1,je
   do i=is+1,ie
    del2(i,j,:)= maskT(i,j,:)* (flux_east(i,j,:) - flux_east(i-1,j,:))/(cost(j)*dxt(i)) &
                              +(flux_north(i,j,:) - flux_north(i,j-1,:))/(cost(j)*dyt(j)) 
  enddo
 enddo
 call border_exchg_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,del2) 
 call setcyclic_xyz   (is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,del2)

 do j=js,je
   fxa = cost(j)**biha_mixing_cosPower
   do i=is,ie-1
    flux_east(i,j,:)=diff*fxa*(del2(i+1,j,:)-del2(i,j,:))/(cost(j)*dxu(i))*maskU(i,j,:)
   enddo
 enddo
 do j=js,je-1
    fxa = cosu(j)**biha_mixing_cosPower
    flux_north(:,j,:)=diff*fxa*(del2(:,j+1,:)-del2(:,j,:))/dyu(j)*maskV(:,j,:)*cosu(j)
 enddo 
 flux_east(ie,:,:)=0.; flux_north(:,je,:)=0.

 ! update tendency
 do k=1,nz
  do j=js_pe,je_pe
   do i=is_pe,ie_pe
      fxa= maskT(i,j,k)*( (flux_east(i,j,k) - flux_east(i-1,j,k))/(cost(j)*dxt(i)) &
                         +(flux_north(i,j,k) - flux_north(i,j-1,k))/(cost(j)*dyt(j))  )
      tra(i,j,k,taup1)=tra(i,j,k,taup1)+dt_tracer*fxa
   enddo
  enddo
 enddo

end subroutine tracer_biharmonic




subroutine tracer_diffusion(is_,ie_,js_,je_,nz_,tra)
!---------------------------------------------------------------------------------
! Diffusion of a tracer, taup1 is only updated but not set
!---------------------------------------------------------------------------------
 use main_module   
 implicit none
 integer :: is_,ie_,js_,je_,nz_
 real*8, dimension(is_:ie_,js_:je_,nz_,3) :: tra
 integer :: i,j,k

  do j=js_pe-onx,je_pe+onx
   do i=is_pe-onx,ie_pe+onx-1
    flux_east(i,j,:)=K_h*(tra(i+1,j,:,tau)-tra(i,j,:,tau))/(cost(j)*dxu(i))*maskU(i,j,:)
   enddo
  enddo
  flux_east(ie_pe+onx,:,:)=0.
  do j=js_pe-onx,je_pe+onx-1
    flux_north(:,j,:)=K_h*(tra(:,j+1,:,tau)-tra(:,j,:,tau))/dyu(j)*maskV(:,j,:)*cosu(j)
  enddo
  flux_north(:,je_pe+onx,:)=0.

  if (enable_hor_friction_cos_scaling) then
   do j=js_pe-onx,je_pe+onx
    flux_east(:,j,:)  = flux_east(:,j,:)*cost(j)**hor_friction_cosPower
    flux_north(:,j,:) = flux_north(:,j,:)*cosu(j)**hor_friction_cosPower
   enddo
  endif

  do k=1,nz
   do j=js_pe,je_pe
    do i=is_pe,ie_pe
       tra(i,j,k,taup1)=tra(i,j,k,taup1)+dt_tracer*maskT(i,j,k)*( &
                      (flux_east(i,j,k) - flux_east(i-1,j,k))/(cost(j)*dxt(i)) &
                     +(flux_north(i,j,k) -flux_north(i,j-1,k))/(cost(j)*dyt(j)) )
    enddo
   enddo
  enddo
end subroutine tracer_diffusion







subroutine tempsalt_sources
!---------------------------------------------------------------------------------
! Sources of temp and salinity, 
! effect on dyn. Enthalpy is stored
!---------------------------------------------------------------------------------
 use main_module   
 implicit none
 integer :: i,j,k,ks
 real*8 :: aloc(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz)

  temp(:,:,:,taup1)=temp(:,:,:,taup1)+dt_tracer*temp_source*maskT
  salt(:,:,:,taup1)=salt(:,:,:,taup1)+dt_tracer*salt_source*maskT

  if (enable_conserve_energy) then
    ! diagnose effect on dynamic enthalpy 
    do k=1,nz
     do j=js_pe-onx+1,je_pe+onx-1
      do i=is_pe-onx+1,ie_pe+onx-1
       aloc(i,j,k) =-grav/rho_0*maskT(i,j,k)*( &
                    int_drhodT(i,j,k,tau)*temp_source(i,j,k)+int_drhodS(i,j,k,tau)*salt_source(i,j,k) )
      enddo
     enddo
    end do
    ! dissipation interpolated on W-grid
    do j=js_pe-onx,je_pe+onx
     do i=is_pe-onx,ie_pe+onx
       ks=kbot(i,j)
       if (ks>0) then
        k=ks; P_diss_sources(i,j,k) = 0.5*(aloc(i,j,k)+aloc(i,j,k+1)) + 0.5*aloc(i,j,k)*dzw(max(1,k-1))/dzw(k)
        do k=ks+1,nz-1
         P_diss_sources(i,j,k) =  0.5*(aloc(i,j,k) +aloc(i,j,k+1))
        enddo
        k=nz; P_diss_sources(i,j,k) = aloc(i,j,k)
       endif
     enddo
    enddo
  endif
end subroutine  tempsalt_sources






