

subroutine biharmonic_thickness_backscatter_energy
 use main_module   
 use biharmonic_thickness_module
 implicit none
 integer :: i,j,k
 real*8 :: forc(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx),fxa
 real*8 :: back(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx)
 real*8 :: taper,z1,z2
 taper(z1,z2) = 0.5*(1.+tanh((-z1-z2)/10.)) 

 forc=0d0
 back=0d0
 do k=1,nz
   do j=js_pe,je_pe
    do i=is_pe,ie_pe 
     forc(i,j) = forc(i,j) + (-diss_biha_skew(i,j,k)+K_diss_h(i,j,k))*dzw(k)*maskW(i,j,k)
     back(i,j) = back(i,j) + diss_back(i,j,k)*dzw(k)*maskW(i,j,k)
    enddo
   enddo
 enddo
  
 ! horizontal diffusion of energy
 do j=js_pe,je_pe
  do i=is_pe-1,ie_pe
   flux_east(i,j,nz)=e_back_diff*(e_back(i+1,j,tau)-e_back(i,j,tau))/(cost(j)*dxu(i))*maskU(i,j,nz)
  enddo
 enddo
 do j=js_pe-1,je_pe
  do i=is_pe,ie_pe
   flux_north(i,j,nz)=e_back_diff*(e_back(i,j+1,tau)-e_back(i,j,tau))/dyu(j)*maskV(i,j,nz)*cosu(j)
  enddo 
 enddo

 do j=js_pe,je_pe
  do i=is_pe,ie_pe  
    fxa = (flux_east(i,j,nz) - flux_east(i-1,j,nz))/(cost(j)*dxt(i)) &
         +(flux_north(i,j,nz) -flux_north(i,j-1,nz))/(cost(j)*dyt(j)) 
    e_back(i,j,taup1) = (e_back(i,j,tau) + dt_tracer*( fxa + forc(i,j) - back(i,j) ) ) &
                       /(1+dt_tracer*e_back_damp*e_back(i,j,tau))*maskW(i,j,nz)
  enddo
 enddo  
 call border_exchg_xy(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,e_back(:,:,taup1)) 
 call setcyclic_xy(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,e_back(:,:,taup1))

 do k=1,nz
  do j=js_pe-onx,je_pe+onx
   do i=is_pe-onx,ie_pe+onx    
    A_thk(i,j,k) = - min(A_thk_max, L_back*sqrt( max(0d0,e_back(i,j,tau)) ) )
   enddo
  enddo 
 enddo
 
 if (enable_biharmonic_thickness_backscatter_taper_mld) then
  do k=1,nz
   do j=js_pe-onx,je_pe+onx
    do i=is_pe-onx,ie_pe+onx  
      A_thk(i,j,k) = A_thk(i,j,k) * taper(zw(k),mld(i,j)+mld0_thk)  
    enddo
   enddo 
  enddo
 endif
 
end subroutine biharmonic_thickness_backscatter_energy



subroutine biharmonic_thickness_backscatter_skew(is_,ie_,js_,je_,nz_,drdTS,ddzt,ddxt,ddyt)
 use main_module   
 use biharmonic_thickness_module
 use diagnostics_module
 implicit none
 integer :: is_,ie_,js_,je_,nz_
 real*8 :: drdTS(is_:ie_,js_:je_,nz_,2)
 real*8 :: ddzt(is_:ie_,js_:je_,nz_,2)
 real*8 :: ddxt(is_:ie_,js_:je_,nz_,2)
 real*8 :: ddyt(is_:ie_,js_:je_,nz_,2)
 
 integer :: i,j,k,ip,jp,kr,n
 real*8 :: Ai_ez(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz,0:1,0:1) 
 real*8 :: Ai_nz(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz,0:1,0:1) 
 real*8 :: Ai_bx(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz,0:1,0:1)
 real*8 :: Ai_by(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz,0:1,0:1) 
 real*8 :: diff_loc(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz)
 real*8 :: drodxe,drodze,drodyn,drodzn,drodxb,drodyb,drodzb
 real*8 :: taper,sxe,syn,sxb,syb,clip,z1,z2
 real*8,parameter :: iso_dslope=0.002,iso_slopec=0.01, epsln=1.D-20 
 taper(sxe,z1,z2) = 0.5*(1.+tanh((-abs(sxe)+iso_slopec)/iso_dslope)) *0.5*(1.+tanh((-z1-z2)/10.)) 
 !taper(sxe) = 0.5*(1.+tanh((-abs(sxe)+iso_slopec)/iso_dslope)) 
!-----------------------------------------------------------------------
 ! statement functions for density triads
!-----------------------------------------------------------------------
 drodxe(i,j,k,ip)     = drdTS(i+ip,j,k,1)*ddxt(i,j,k,1)          + drdTS(i+ip,j,k,2)*ddxt(i,j,k,2) 
 drodze(i,j,k,ip,kr)  = drdTS(i+ip,j,k,1)*ddzt(i+ip,j,k+kr-1,1)  + drdTS(i+ip,j,k,2)*ddzt(i+ip,j,k+kr-1,2)
 drodyn(i,j,k,jp)     = drdTS(i,j+jp,k,1)*ddyt(i,j,k,1)          + drdTS(i,j+jp,k,2)*ddyt(i,j,k,2) 
 drodzn(i,j,k,jp,kr)  = drdTS(i,j+jp,k,1)*ddzt(i,j+jp,k+kr-1,1)  + drdTS(i,j+jp,k,2)*ddzt(i,j+jp,k+kr-1,2)
 drodxb(i,j,k,ip,kr)  = drdTS(i,j,k+kr,1)*ddxt(i-1+ip,j,k+kr,1)  + drdTS(i,j,k+kr,2)*ddxt(i-1+ip,j,k+kr,2) 
 drodyb(i,j,k,jp,kr)  = drdTS(i,j,k+kr,1)*ddyt(i,j-1+jp,k+kr,1)  + drdTS(i,j,k+kr,2)*ddyt(i,j-1+jp,k+kr,2) 
 drodzb(i,j,k,kr)     = drdTS(i,j,k+kr,1)*ddzt(i,j,k,1)          + drdTS(i,j,k+kr,2)*ddzt(i,j,k,2)

!-----------------------------------------------------------------------
! smooth the gradients
!-----------------------------------------------------------------------
 do n=1,2
  if (biharmonic_thickness_backscatter_smooth == -1) then
   ! Do nothing
  else if (biharmonic_thickness_backscatter_smooth == 0) then
    call biharmonic_thickness_smooth_1D_1(ddzt(:,:,:,n),maskW,1)
    ddzt(:,:,nz,:)=0
  else if (biharmonic_thickness_backscatter_smooth == 1) then
    call biharmonic_thickness_smooth_1D_2(ddzt(:,:,:,n),maskW,1)
    call biharmonic_thickness_smooth_1D_2(ddxt(:,:,:,n),maskU,1)
    call biharmonic_thickness_smooth_1D_2(ddyt(:,:,:,n),maskV,1)
    ddzt(:,:,nz,:)=0 
  else if (biharmonic_thickness_backscatter_smooth == 3) then
    call border_exchg_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,ddxt(:,:,:,n)) 
    call setcyclic_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,ddxt(:,:,:,n)) 
    call border_exchg_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,ddyt(:,:,:,n)) 
    call setcyclic_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,ddyt(:,:,:,n))  
    !call biharmonic_thickness_smooth(drdTS(:,:,:,n),maskT,1)
    call biharmonic_thickness_smooth_3D_2(ddzt(:,:,:,n),maskW,1)
    call biharmonic_thickness_smooth_3D_2(ddxt(:,:,:,n),maskU,1)
    call biharmonic_thickness_smooth_3D_2(ddyt(:,:,:,n),maskV,1)
    ddzt(:,:,nz,:)=0
  else
   call halt_stop('ERROR (1) in biharmonic_thickness_backscatter_skew')
  endif 
 enddo


 if (enable_biharmonic_thickness_backscatter_integrate_energy) then
  call biharmonic_thickness_backscatter_energy
 else
  A_thk=A_thk_0
 endif

 Ai_ez=0d0;Ai_nz=0d0;Ai_bx=0d0;Ai_by=0d0
!-----------------------------------------------------------------------
!     Compute Ai_ez on center of east face of T cell.
!-----------------------------------------------------------------------
 do k=2,nz
  do j=js_pe,je_pe
   do i=is_pe-1,ie_pe
    do kr=0,1
     do ip=0,1
       sxe  = -drodxe(i,j,k,ip)/(min(0d0,drodze(i,j,k,ip,kr))-epsln)  ! i+1, k-1
       clip = taper(sxe,zw(k),mldu(i,j)+mld0_thk)
       Ai_ez(i,j,k,ip,kr) =  clip*sxe*maskU(i,j,k)
     enddo
    enddo
   enddo
  enddo
 enddo
 k=1
 do j=js_pe,je_pe
   do i=is_pe-1,ie_pe
    kr=1
    do ip=0,1
       sxe  = -drodxe(i,j,k,ip)/(min(0d0,drodze(i,j,k,ip,kr))-epsln)
       clip = taper(sxe,zw(k),mldu(i,j)+mld0_thk)
       Ai_ez(i,j,k,ip,kr) =  clip*sxe*maskU(i,j,k)
    enddo
   enddo
 enddo
 Ai_ez(:,:,1,:,0) =0d0
!-----------------------------------------------------------------------
!     Compute Ai_nz on center of north face of T cell.
!-----------------------------------------------------------------------
 do k=2,nz
  do j=js_pe-1,je_pe
   do i=is_pe,ie_pe
    do kr=0,1
      do jp=0,1
         syn = -drodyn(i,j,k,jp)/(min(0d0,drodzn(i,j,k,jp,kr))-epsln)
         clip = taper(syn,zw(k),mldv(i,j)+mld0_thk)
         Ai_nz(i,j,k,jp,kr) = clip*syn*maskV(i,j,k)
      enddo
     enddo
    enddo
  enddo
 enddo
 k=1
 do j=js_pe-1,je_pe
   do i=is_pe,ie_pe
     kr=1
     do jp=0,1
         syn = -drodyn(i,j,k,jp)/(min(0d0,drodzn(i,j,k,jp,kr))-epsln)
         clip = taper(syn,zw(k),mldv(i,j)+mld0_thk)
         Ai_nz(i,j,k,jp,kr) = clip*syn*maskV(i,j,k)         
     enddo
   enddo
 enddo
 Ai_nz(:,:,1,:,0) = 0d0 
!-----------------------------------------------------------------------
! compute Ai_bx, Ai_by on top face of T cell.
!-----------------------------------------------------------------------
 do k=1,nz-1
  do j=js_pe,je_pe
   do i=is_pe,ie_pe
!    eastward slopes at the top of T cells
     do ip=0,1
      do kr=0,1 
        sxb = -drodxb(i,j,k,ip,kr)/(min(0d0,drodzb(i,j,k,kr))-epsln)  ! i-1,k+1
        clip = taper(sxb,zw(k),mld(i,j)+mld0_thk)
        Ai_bx(i,j,k,ip,kr) =  clip*sxb*maskW(i,j,k)
      enddo
     enddo
!    northward slopes at the top of T cells
     do jp=0,1
      do kr=0,1
        syb = -drodyb(i,j,k,jp,kr)/(min(0d0,drodzb(i,j,k,kr))-epsln)
        clip = taper(syb,zw(k),mld(i,j)+mld0_thk)
        Ai_by(i,j,k,jp,kr) = clip*syb  *maskW(i,j,k)
      enddo
     enddo
   enddo
  enddo
 enddo
 Ai_bx(:,:,nz,:,:)=0d0
 Ai_by(:,:,nz,:,:)=0d0
 

 if (enable_conserve_energy) diss_back = 0d0
 
 do j=js_pe-onx,je_pe+onx
  diff_loc(:,j,:) = -A_thk(:,j,:) ! since A_thkbi>0
 enddo 
 
 call biharmonic_thickness_skew_divergence(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,&
                                           temp,Ai_ez,Ai_nz,Ai_bx,Ai_by,diff_loc,diss_back, .true.)
 call biharmonic_thickness_skew_divergence(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz, &
                                           salt,Ai_ez,Ai_nz,Ai_bx,Ai_by,diff_loc,diss_back, .false.)

 if (mod(itt,int(snapint/dt_tracer)) ==0 .or. itt==0) &
             call biharmonic_thickness_diag_store_tensor_back(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz, &
                                                         -Ai_by(:,:,:,0,0)*diff_loc,Ai_bx(:,:,:,0,0)*diff_loc)

                                                                                           
end subroutine biharmonic_thickness_backscatter_skew





subroutine biharmonic_thickness_backscatter_harmonic
 use main_module
 use biharmonic_thickness_module
 implicit none
 integer :: i,j,k
 real*8 :: diss(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz),fxa

 if (enable_biharmonic_thickness_backscatter_integrate_energy) then
  call biharmonic_thickness_backscatter_energy
 else
  A_thk=A_thk_0
 endif

 do k=1,nz
  do j=js_pe-onx,je_pe+onx
   do i=is_pe-onx,ie_pe+onx-1
    fxa = A_thk(i+1,j,k)
    flux_east(i,j,k) = fxa*(u(i+1,j,k,tau)-u(i,j,k,tau))/(cost(j)*dxt(i+1))*maskU(i+1,j,k)*maskU(i,j,k)
   enddo
  enddo
  do j=js_pe-onx,je_pe+onx-1
   do i=is_pe-onx,ie_pe+onx-1   
    fxa = 0.25*(A_thk(i,j,k) + A_thk(i,j+1,k) + A_thk(i+1,j,k) + A_thk(i+1,j+1,k))
    flux_north(i,j,k) = fxa*(u(i,j+1,k,tau)-u(i,j,k,tau))/dyu(j)*maskU(i,j+1,k)*maskU(i,j,k)*cosu(j)
   enddo 
  enddo 
 enddo
 

 do j=js_pe,je_pe
   do i=is_pe,ie_pe
    du_mix(i,j,:)= du_mix(i,j,:) + maskU(i,j,:)*((flux_east(i,j,:) - flux_east(i-1,j,:))/(cost(j)*dxu(i)) &
                                                +(flux_north(i,j,:) - flux_north(i,j-1,:))/(cost(j)*dyt(j)) )
  enddo
 enddo

 if (enable_conserve_energy) then
  diss_back = 0d0
  do k=1,nz
   do j=js_pe,je_pe
    do i=is_pe-1,ie_pe
     diss(i,j,k) =-0.5*((u(i+1,j,k,tau)-u(i,j,k,tau))*flux_east(i,j,k) &
                       +(u(i,j,k,tau)-u(i-1,j,k,tau))*flux_east(i-1,j,k))/(cost(j)*dxu(i))  &
                  -0.5*((u(i,j+1,k,tau)-u(i,j,k,tau))*flux_north(i,j,k)+ &
                        (u(i,j,k,tau)-u(i,j-1,k,tau))*flux_north(i,j-1,k))/(cost(j)*dyt(j)) 
    enddo
   enddo
  enddo
  call calc_diss(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,diss,diss_back,'U')
 endif

 do k=1,nz
  do j=js_pe-onx,je_pe+onx-1
   do i=is_pe-onx,ie_pe+onx-1    
    fxa = 0.25*(A_thk(i,j,k) + A_thk(i,j+1,k) + A_thk(i+1,j,k) + A_thk(i+1,j+1,k))
    flux_east(i,j,k) = fxa*(v(i+1,j,k,tau)-v(i,j,k,tau))/(cosu(j)*dxu(i)) *maskV(i+1,j,k)*maskV(i,j,k)
   enddo
  enddo
  do j=js_pe-onx,je_pe+onx-1
   do i=is_pe-onx,ie_pe+onx
    fxa = A_thk(i,j+1,k)
    flux_north(i,j,k) = fxa*(v(i,j+1,k,tau)-v(i,j,k,tau) )/dyt(j+1)*cost(j+1)*maskV(i,j,k)*maskV(i,j+1,k)
   enddo 
  enddo
 enddo
 
 do j=js_pe,je_pe
   do i=is_pe,ie_pe
    dv_mix(i,j,:)= dv_mix(i,j,:) + maskV(i,j,:)*( (flux_east(i,j,:) - flux_east(i-1,j,:))/(cosu(j)*dxt(i))  &
                                                 +(flux_north(i,j,:) - flux_north(i,j-1,:))/(dyu(j)*cosu(j)) )
  enddo
 enddo

 if (enable_conserve_energy) then
  do k=1,nz
   do j=js_pe-1,je_pe
    do i=is_pe,ie_pe
     diss(i,j,k) =-0.5*((v(i+1,j,k,tau)-v(i,j,k,tau))*flux_east(i,j,k)+ &
                        (v(i,j,k,tau)-v(i-1,j,k,tau))*flux_east(i-1,j,k))/(cosu(j)*dxt(i)) &
                  -0.5*((v(i,j+1,k,tau)-v(i,j,k,tau))*flux_north(i,j,k)+ &
                        (v(i,j,k,tau)-v(i,j-1,k,tau))*flux_north(i,j-1,k))/(cosu(j)*dyu(j)) 
    enddo
   enddo
  enddo
  call calc_diss(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,diss,diss_back,'V')
 endif
end subroutine biharmonic_thickness_backscatter_harmonic






subroutine biharmonic_thickness_smooth_3D_2(var,mask_loc,m)
 use main_module   
 implicit none
 integer :: i,j,k,ii,jj,n,m,kk,kkk
 real*8 :: var(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz)
 real*8 :: mask_loc(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz)
 real*8 :: aloc(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz)
 real*8 :: bloc(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz)
 
 do n=1,m
  aloc = 0d0
  bloc = 0d0
  do k=1,nz
do kk=-2,2
 kkk = min(nz,max(1,k+kk))
   do ii=-2,2
    do jj=-2,2
     do j=js_pe,je_pe
      do i=is_pe,ie_pe 
       aloc(i,j,k) =  aloc(i,j,k) + var(i+ii,j+jj,kkk)*mask_loc(i+ii,j+jj,kkk) 
       bloc(i,j,k) =  bloc(i,j,k) + mask_loc(i+ii,j+jj,kkk) 
      enddo
     enddo     
    enddo
   enddo 
  enddo 
enddo  
  where (bloc/=0d0) var = aloc/bloc
  
  call border_exchg_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,var) 
  call setcyclic_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,var) 
 enddo 
end subroutine biharmonic_thickness_smooth_3D_2



subroutine biharmonic_thickness_smooth_3D_1(var,mask_loc,m)
 use main_module   
 implicit none
 integer :: i,j,k,ii,jj,n,m,kk,kkk
 real*8 :: var(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz)
 real*8 :: mask_loc(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz)
 real*8 :: aloc(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz)
 real*8 :: bloc(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz)
 
 do n=1,m
  aloc = 0d0
  bloc = 0d0
  do k=1,nz
do kk=-1,1
 kkk = min(nz,max(1,k+kk))
   do ii=-1,1
    do jj=-1,1
     do j=js_pe,je_pe
      do i=is_pe,ie_pe 
       aloc(i,j,k) =  aloc(i,j,k) + var(i+ii,j+jj,kkk)*mask_loc(i+ii,j+jj,kkk) 
       bloc(i,j,k) =  bloc(i,j,k) + mask_loc(i+ii,j+jj,kkk) 
      enddo
     enddo     
    enddo
   enddo 
  enddo 
enddo  
  where (bloc/=0d0) var = aloc/bloc
  
  call border_exchg_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,var) 
  call setcyclic_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,var) 
 enddo 
end subroutine biharmonic_thickness_smooth_3D_1




subroutine biharmonic_thickness_smooth_1D_1(var,mask_loc,m)
 use main_module   
 implicit none
 integer :: i,j,k,n,m,kk,kkk
 real*8 :: var(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz)
 real*8 :: mask_loc(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz)
 real*8 :: aloc(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz)
 real*8 :: bloc(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz)
 
 do n=1,m
  aloc = 0d0
  bloc = 0d0
  do k=1,nz
do kk=-1,1
 kkk = min(nz,max(1,k+kk))
     do j=js_pe,je_pe
      do i=is_pe,ie_pe 
       aloc(i,j,k) =  aloc(i,j,k) + var(i,j,kkk)*mask_loc(i,j,kkk) 
       bloc(i,j,k) =  bloc(i,j,k) + mask_loc(i,j,kkk) 
      enddo
     enddo 
  enddo 
enddo  
  where (bloc/=0d0) var = aloc/bloc
  enddo 
end subroutine biharmonic_thickness_smooth_1D_1


subroutine biharmonic_thickness_smooth_1D_2(var,mask_loc,m)
 use main_module   
 implicit none
 integer :: i,j,k,n,m,kk,kkk
 real*8 :: var(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz)
 real*8 :: mask_loc(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz)
 real*8 :: aloc(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz)
 real*8 :: bloc(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz)
 
 do n=1,m
  aloc = 0d0
  bloc = 0d0
  do k=1,nz
do kk=-2,2
 kkk = min(nz,max(1,k+kk))
     do j=js_pe,je_pe
      do i=is_pe,ie_pe 
       aloc(i,j,k) =  aloc(i,j,k) + var(i,j,kkk)*mask_loc(i,j,kkk) 
       bloc(i,j,k) =  bloc(i,j,k) + mask_loc(i,j,kkk) 
      enddo
     enddo 
  enddo 
enddo  
  where (bloc/=0d0) var = aloc/bloc
  enddo 
end subroutine biharmonic_thickness_smooth_1D_2





subroutine biharmonic_thickness_smooth_2D_2(var,mask_loc,m)
 use main_module   
 implicit none
 integer :: i,j,k,ii,jj,n,m
 real*8 :: var(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz)
 real*8 :: mask_loc(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz)
 real*8 :: aloc(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz)
 real*8 :: bloc(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz)
 
 do n=1,m
  aloc = 0d0
  bloc = 0d0
  do k=1,nz

   do ii=-onx,onx
    do jj=-onx,onx
     do j=js_pe,je_pe
      do i=is_pe,ie_pe 
       aloc(i,j,k) =  aloc(i,j,k) + var(i+ii,j+jj,k)*mask_loc(i+ii,j+jj,k) 
       bloc(i,j,k) =  bloc(i,j,k) + mask_loc(i+ii,j+jj,k) 
      enddo
     enddo     
    enddo
   enddo 
  enddo 

  where (bloc/=0d0) var = aloc/bloc
  
  call border_exchg_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,var) 
  call setcyclic_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,var) 
 enddo 
end subroutine biharmonic_thickness_smooth_2D_2




subroutine biharmonic_thickness_smooth_2D_1(var,mask_loc,m)
 use main_module   
 implicit none
 integer :: i,j,k,ii,jj,n,m
 real*8 :: var(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz)
 real*8 :: mask_loc(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz)
 real*8 :: aloc(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz)
 real*8 :: bloc(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz)
 
 do n=1,m
  aloc = 0d0
  bloc = 0d0
  do k=1,nz

   do ii=-1,1
    do jj=-1,1
     do j=js_pe,je_pe
      do i=is_pe,ie_pe 
       aloc(i,j,k) =  aloc(i,j,k) + var(i+ii,j+jj,k)*mask_loc(i+ii,j+jj,k) 
       bloc(i,j,k) =  bloc(i,j,k) + mask_loc(i+ii,j+jj,k) 
      enddo
     enddo     
    enddo
   enddo 
  enddo 

  where (bloc/=0d0) var = aloc/bloc
  
  call border_exchg_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,var) 
  call setcyclic_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,var) 
 enddo 
end subroutine biharmonic_thickness_smooth_2D_1



