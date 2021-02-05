




subroutine rossmix_eddy_streamfunction
 !---------------------------------------------------------------------------------
 ! streamfunction for eddy driven flow 
 !---------------------------------------------------------------------------------
 use main_module   
 use rossmix_module   
 implicit none
 integer :: n,m,k
 real*8 :: Azloc(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz)
 real*8 :: aloc(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx)
 real*8 :: bloc(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx)

 Bx=0; By=0
 do m=1,nmodes

  aloc = 0; bloc = 0
  do n=2,nphi-1
   aloc(:,:) = aloc(:,:) + dphit(n)*G_s(:,:,n,tau,m)*cos(phit(n))*maskTp(:,:,n)
   bloc(:,:) = bloc(:,:) + dphit(n)*G_s(:,:,n,tau,m)*sin(phit(n))*maskTp(:,:,n)
  enddo
  do k=1,nz
   Azloc(:,:,k) = 2*phi_ga_s(:,:,k,m)/sqrt(max(rossmix_N2min,Nsqr(:,:,k,tau)))/max(1d-12,Rn(:,:,m))
  enddo
  do k=1,nz-1
    Bx(:,:,k) = Bx(:,:,k) - Azloc(:,:,k)*aloc(:,:)
    By(:,:,k) = By(:,:,k) - Azloc(:,:,k)*bloc(:,:)
  enddo

  aloc = 0; bloc = 0
  do n=2,nphi-1
   aloc(:,:) = aloc(:,:) + dphit(n)*G_l(:,:,n,tau,m)*cos(phit(n))*maskTp(:,:,n)
   bloc(:,:) = bloc(:,:) + dphit(n)*G_l(:,:,n,tau,m)*sin(phit(n))*maskTp(:,:,n)
  enddo
  do k=1,nz
   Azloc(:,:,k) = 2*phi_ga_l(:,:,k,m)/sqrt(max(rossmix_N2min,Nsqr(:,:,k,tau)))/max(1d-12,Rn(:,:,m))
  enddo
  do k=1,nz-1
    Bx(:,:,k) = Bx(:,:,k) - Azloc(:,:,k)*aloc(:,:)
    By(:,:,k) = By(:,:,k) - Azloc(:,:,k)*bloc(:,:)
  enddo

 enddo

 call tgrid_to_ugrid( is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,By,By)  
 call tgrid_to_vgrid( is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,Bx,Bx)  

end subroutine rossmix_eddy_streamfunction






subroutine rossmix_eddy_advect
 use main_module   
 use rossmix_module   
 implicit none
 integer :: i,j
 real*8 :: dtr(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz)

 call rossmix_eddy_advect_tr(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,temp(:,:,:,tau),dtr)
 do j=js_pe,je_pe
   do i=is_pe,ie_pe
     temp(i,j,:,taup1)=temp(i,j,:,taup1)+dt_tracer*dtr(i,j,:)
   enddo
 enddo

 call rossmix_eddy_advect_tr(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,salt(:,:,:,tau),dtr)
 do j=js_pe,je_pe
   do i=is_pe,ie_pe
     salt(i,j,:,taup1)=salt(i,j,:,taup1)+dt_tracer*dtr(i,j,:)
   enddo
 enddo

 ! calculate dissipation by eddy advection -> in rossmix_0.1

end subroutine rossmix_eddy_advect




subroutine rossmix_eddy_advect_tr(is_,ie_,js_,je_,nz_,var,dvar)
 use main_module   
 use rossmix_module   
 implicit none
 integer :: is_,ie_,js_,je_,nz_
 real*8 :: var(is_:ie_,js_:je_,nz_),dvar(is_:ie_,js_:je_,nz_)
 integer :: i,j,k

 do i=is_pe-1,ie_pe
   flux_east(i,:,:) = 0.5*(var(i,:,:)+var(i+1,:,:))*rossmix_ue(i,:,:)*maskU(i,:,:)
 enddo
 do j=js_pe-1,je_pe
   flux_north(:,j,:) = cosu(j)*0.5*(var(:,j,:)+var(:,j+1,:))*rossmix_ve(:,j,:)*maskV(:,j,:)
 enddo
 do k=1,nz-1
   flux_top(:,:,k)=0.5*(var(:,:,k)+var(:,:,k+1))*rossmix_we(:,:,k)*maskW(:,:,k)
 enddo
 flux_top(:,:,nz)=0.0
 
 do j=js_pe,je_pe
   do i=is_pe,ie_pe
      dvar(i,j,:)=maskT(i,j,:)* (-( flux_east(i,j,:)-  flux_east(i-1,j,:))/(cost(j)*dxt(i)) &
                                  -(flux_north(i,j,:)- flux_north(i,j-1,:))/(cost(j)*dyt(j)) )
   enddo
 enddo
 k=1; dvar(:,:,k)=dvar(:,:,k)-maskT(:,:,k)*flux_top(:,:,k)/dzt(k)   
 do k=2,nz
   dvar(:,:,k)=dvar(:,:,k)-maskT(:,:,k)*(flux_top(:,:,k)- flux_top(:,:,k-1))/dzt(k)   
 enddo
end subroutine rossmix_eddy_advect_tr



subroutine rossmix_eddy_velocity
 use main_module   
 use rossmix_module   
 implicit none
 integer :: i,j,k

 k=1
 rossmix_ue(:,:,k) = -(By(:,:,k))/dzt(k)
 rossmix_ve(:,:,k) =  (Bx(:,:,k))/dzt(k)
 do k=2,nz
   rossmix_ue(:,:,k) = -(By(:,:,k)-By(:,:,k-1))/dzt(k)
   rossmix_ve(:,:,k) =  (Bx(:,:,k)-Bx(:,:,k-1))/dzt(k)
 enddo
 rossmix_ue = rossmix_ue*maskU
 rossmix_ve = rossmix_ve*maskV

 call border_exchg_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,rossmix_ue) 
 call setcyclic_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,rossmix_ue)
 call border_exchg_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,rossmix_ve) 
 call setcyclic_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,rossmix_ve)

 k=1
 do j=js_pe-onx+1,je_pe+onx
   do i=is_pe-onx+1,ie_pe+onx
     rossmix_we(i,j,k) =-maskW(i,j,k)*dzt(k)* &
               ((        rossmix_ue(i,j,k)-          rossmix_ue(i-1,j,k))/(cost(j)*dxt(i)) &
               +(cosu(j)*rossmix_ve(i,j,k)-cosu(j-1)*rossmix_ve(i,j-1,k))/(cost(j)*dyt(j)) )
   enddo
 enddo
 do k=2,nz
  do j=js_pe-onx+1,je_pe+onx
   do i=is_pe-onx+1,ie_pe+onx
     rossmix_we(i,j,k) = rossmix_we(i,j,k-1)-maskW(i,j,k)*dzt(k)* &
               ((        rossmix_ue(i,j,k)          -rossmix_ue(i-1,j,k))/(cost(j)*dxt(i)) &
               +(cosu(j)*rossmix_ve(i,j,k)-cosu(j-1)*rossmix_ve(i,j-1,k))/(cost(j)*dyt(j)) )
   enddo
  enddo
 enddo
end subroutine rossmix_eddy_velocity







subroutine rossmix_lateral_stress
 use main_module   
 use rossmix_module   
 implicit none
 integer :: i,j,k,m,n
 real*8 :: Cloc(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx)
 real*8 :: Dloc(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx)

 !---------------------------------------------------------------------------------
 ! lateral eddy zonal and meridional momentum fluxes
 !---------------------------------------------------------------------------------
 stress_ux = 0; stress_uy = 0; stress_vx = 0; stress_vy = 0
 do m=1,nmodes
   Cloc = 0.0; Dloc = 0.0
   do n=2,nphi-1
     Cloc(:,:) = Cloc(:,:) + 2*dphit(n)*(E_s(:,:,n,tau,m)*phi_nu_s(:,:,m)+E_l(:,:,n,tau,m)*phi_nu_l(:,:,m)) &
                                                    *maskTp(:,:,n)*sin(phit(n))**2
     Dloc(:,:) = Dloc(:,:) - 2*dphit(n)*(E_s(:,:,n,tau,m)*phi_nu_s(:,:,m)+E_l(:,:,n,tau,m)*phi_nu_l(:,:,m)) &
                                                    *maskTp(:,:,n)*sin(phit(n))*cos(phit(n))
   enddo
   do k=1,nz
    stress_ux(:,:,k) = stress_ux(:,:,k) + Cloc(:,:)*maskT(:,:,k)*phin(:,:,k,m)**2
    stress_uy(:,:,k) = stress_uy(:,:,k) + Dloc(:,:)*maskT(:,:,k)*phin(:,:,k,m)**2
   enddo

   Cloc = 0.0;
   do n=2,nphi-1
     Cloc(:,:) = Cloc(:,:) + 2*dphit(n)*(E_s(:,:,n,tau,m)*phi_nu_s(:,:,m)+E_l(:,:,n,tau,m)*phi_nu_l(:,:,m)) &
                                                    *maskTp(:,:,n)*cos(phit(n))**2
   enddo
   do k=1,nz
    stress_vx(:,:,k) = stress_vx(:,:,k) + Dloc(:,:)*maskT(:,:,k)*phin(:,:,k,m)**2
    stress_vy(:,:,k) = stress_vy(:,:,k) + Cloc(:,:)*maskT(:,:,k)*phin(:,:,k,m)**2
   enddo
 enddo

 !---------------------------------------------------------------------------------
 ! Interpolation
 !---------------------------------------------------------------------------------
 stress_ux(is_pe-onx:ie_pe+onx-1,:,:) = stress_ux(is_pe-onx+1:ie_pe+onx,:,:)
 call tgrid_to_ugrid( is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,stress_uy,stress_uy)  
 call tgrid_to_vgrid( is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,stress_uy,stress_uy)  

 stress_vy(is_pe-onx:ie_pe+onx-1,:,:) = stress_vy(is_pe-onx+1:ie_pe+onx,:,:)
 call tgrid_to_ugrid( is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,stress_vx,stress_vx)  
 call tgrid_to_vgrid( is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,stress_vx,stress_vx)  

 !---------------------------------------------------------------------------------
 ! Masking
 !---------------------------------------------------------------------------------
 do i=is_pe-onx,ie_pe+onx-1
   stress_ux(i,:,:) = stress_ux(i,:,:)*maskU(i,:,:)*maskU(i+1,:,:)
 enddo
 do j=js_pe-onx,je_pe+onx-1
    stress_uy(:,j,:) = stress_uy(:,j,:)*maskU(:,j,:)*maskU(:,j+1,:)
 enddo
 do i=is_pe-onx,ie_pe+onx-1
   stress_vx(i,:,:) = stress_vx(i,:,:)*maskV(i,:,:)*maskV(i+1,:,:)
 enddo
 do j=js_pe-onx,je_pe+onx-1
   stress_vy(:,j,:) = stress_vy(:,j,:)*maskV(:,j,:)*maskV(:,j+1,:)
 enddo

end subroutine rossmix_lateral_stress





subroutine rossmix_friction
!=======================================================================
!  rossby wave drag 
!=======================================================================
 use main_module   
 use rossmix_module   
 implicit none
 integer :: i,j,k
 real*8 :: diss(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz)
 !real*8 :: mdiss,fxa


 if (enable_rossmix_lateral_stress) then 

 K_diss_ross = 0.
 !---------------------------------------------------------------------------------
 ! divergence of horizontal eddy zonal momentum flux 
 !---------------------------------------------------------------------------------
 flux_east(:,:,:) = -stress_ux(:,:,:)
 do j=js_pe-1,je_pe
   flux_north(:,j,:)= -stress_uy(:,j,:)*cosu(j)
 enddo
 do j=js_pe,je_pe
   do i=is_pe,ie_pe
    du_mix(i,j,:)= du_mix(i,j,:) + maskU(i,j,:)*((flux_east(i,j,:) -   flux_east(i-1,j,:))/(cost(j)*dxu(i)) &
                                                +(flux_north(i,j,:) - flux_north(i,j-1,:))/(cost(j)*dyt(j)) )
  enddo
 enddo
 !---------------------------------------------------------------------------------
 ! diagnose dissipation by that flux
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
  call calc_diss(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,diss,K_diss_ross,'U')

 !---------------------------------------------------------------------------------
 ! divergence of horizontal eddy meridional momentum flux 
 !---------------------------------------------------------------------------------
 flux_east(:,:,:) = -stress_vx(:,:,:)
 do j=js_pe-1,je_pe
   flux_north(:,j,:)= -stress_vy(:,j,:)*cost(j+1)
 enddo
 do j=js_pe,je_pe
   do i=is_pe,ie_pe
    dv_mix(i,j,:)= dv_mix(i,j,:) + maskV(i,j,:)*( (flux_east(i,j,:) -   flux_east(i-1,j,:))/(cosu(j)*dxt(i))  &
                                                 +(flux_north(i,j,:) - flux_north(i,j-1,:))/(dyu(j)*cosu(j)) )
  enddo
 enddo

 !---------------------------------------------------------------------------------
 ! diagnose dissipation by that flux
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
  call calc_diss(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,diss,K_diss_ross,'V')
 endif

end subroutine rossmix_friction







