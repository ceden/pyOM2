


subroutine integrate_rossmix
 use main_module   
 use rossmix_module   
 use timing_module   
 implicit none
 integer :: i,j,n,m,l,nmax,tauloc
 real*8 :: dtloc,fxa
 real*8 :: advfe_s(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nphi) 
 real*8 :: advfn_s(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nphi) 
 real*8 :: advft_s(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nphi) 
 real*8 :: advfe_l(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nphi) 
 real*8 :: advfn_l(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nphi) 
 real*8 :: advft_l(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nphi) 

 real*8 :: growthE_l(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nphi) 
 real*8 :: growthE_s(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nphi) 
 real*8 :: aloc(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nphi) 

 do m=1,nmodes

  !---------------------------------------------------------------------------------
  ! Diffusion in x,y, and phi splitted from advection
  !---------------------------------------------------------------------------------

  if (enable_rossmix_diffusion) then
   call diff_flux_rossmix(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nphi,advfe_s,advfn_s,advft_s,E_s(:,:,:,:,m),tau)
   call diff_flux_rossmix(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nphi,advfe_l,advfn_l,advft_l,E_l(:,:,:,:,m),tau)
   do n=2,nphi-1
    do j=js_pe,je_pe
     do i=is_pe,ie_pe
      fxa=maskTp(i,j,n)*((advfe_s(i,j,n)- advfe_s(i-1,j,n))/(cost(j)*dxt(i)) &
                        +(advfn_s(i,j,n)- advfn_s(i,j-1,n))/(cost(j)*dyt(j)) &
                        +(advft_s(i,j,n)- advft_s(i,j,n-1))/dphit(n)         )    
      E_s(i,j,n,taup1,m)=E_s(i,j,n,tau,m)  +dt_tracer*fxa 

      fxa=maskTp(i,j,n)*((advfe_l(i,j,n)- advfe_l(i-1,j,n))/(cost(j)*dxt(i)) &
                        +(advfn_l(i,j,n)- advfn_l(i,j-1,n))/(cost(j)*dyt(j)) &
                        +(advft_l(i,j,n)- advft_l(i,j,n-1))/dphit(n)         )    
      E_l(i,j,n,taup1,m)=E_l(i,j,n,tau,m) + dt_tracer*fxa 
     enddo
    enddo
   enddo
  else
   E_l(:,:,:,taup1,m) = E_l(:,:,:,tau,m) 
   E_s(:,:,:,taup1,m) = E_s(:,:,:,tau,m) 
  endif

  !---------------------------------------------------------------------------------
  ! Advection in x,y, and phi and linear damping
  !---------------------------------------------------------------------------------

  tauloc = tau
  dtloc = dt_tracer; nmax=1; 
  if (m==1) then ! time step splitting for barotropic mode
     dtloc = dt_tracer/rossmix_barotropic_fac ; nmax=rossmix_barotropic_fac
  endif

  !---------------------------------------------------------------------------------
  ! net linear growth rate
  !---------------------------------------------------------------------------------

  growthE_l(:,:,:) = lambda_l(:,:,:,m) + 2*nu_l(:,:,:,m)
  growthE_s(:,:,:) = lambda_s(:,:,:,m) + 2*nu_s(:,:,:,m)
  if (m==1) then
   do n=2,nphi-1
    growthE_l(:,:,n) = growthE_l(:,:,n)  - ray_fric_l(:,:,m) 
    growthE_s(:,:,n) = growthE_s(:,:,n)  - ray_fric_s(:,:,m) 
   enddo
  endif

  !---------------------------------------------------------------------------------
  ! loop over time splitting index
  !---------------------------------------------------------------------------------
  do l=1,nmax

   !---------------------------------------------------------------------------------
   ! advective fluxes
   !---------------------------------------------------------------------------------
   call adv_flux_rossmix(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nphi,advfe_s,advfn_s,advft_s,E_s(:,:,:,:,m),tauloc,dtloc, &
                          cgu_s(:,:,:,m),cgv_s(:,:,:,m),phidot_s(:,:,:,m))
   call adv_flux_rossmix(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nphi,advfe_l,advfn_l,advft_l,E_l(:,:,:,:,m),tauloc,dtloc, &
                          cgu_l(:,:,:,m),cgv_l(:,:,:,m),phidot_l(:,:,:,m))
   !---------------------------------------------------------------------------------
   ! reflection of fluxes
   !---------------------------------------------------------------------------------
   call reflect_rossmix(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nphi,m,advfe_s,advfn_s,advfe_l,advfn_l, &
                         E_l(:,:,:,tauloc,m),E_s(:,:,:,tauloc,m))

   !---------------------------------------------------------------------------------
   ! divergence of fluxes and linear damping terms
   !---------------------------------------------------------------------------------
   do n=2,nphi-1
    do j=js_pe,je_pe
     do i=is_pe,ie_pe
      fxa=maskTp(i,j,n)*(-(advfe_s(i,j,n)- advfe_s(i-1,j,n))/(cost(j)*dxt(i)) &
                         -(advfn_s(i,j,n)- advfn_s(i,j-1,n))/(cost(j)*dyt(j)) &
                         -(advft_s(i,j,n)- advft_s(i,j,n-1))/dphit(n)         )    
      E_s(i,j,n,taup1,m)=E_s(i,j,n,taup1,m)  +dtloc*fxa  + dtloc*growthE_s(i,j,n)*E_s(i,j,n,tauloc,m)*maskTp(i,j,n)

      fxa=maskTp(i,j,n)*(-(advfe_l(i,j,n)- advfe_l(i-1,j,n))/(cost(j)*dxt(i)) &
                         -(advfn_l(i,j,n)- advfn_l(i,j-1,n))/(cost(j)*dyt(j)) &
                         -(advft_l(i,j,n)- advft_l(i,j,n-1))/dphit(n)         )    
      E_l(i,j,n,taup1,m)=E_l(i,j,n,taup1,m) + dtloc*fxa  + dtloc*growthE_l(i,j,n)*E_l(i,j,n,tauloc,m)*maskTp(i,j,n)
     enddo
    enddo
   enddo

   !---------------------------------------------------------------------------------
   ! boundary exchange over PEs and cyclic boundary condition
   !---------------------------------------------------------------------------------
   call border_exchg_xyp(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nphi,E_s(:,:,:,taup1,m)) 
   call setcyclic_xyp(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nphi,E_s(:,:,:,taup1,m))
   call border_exchg_xyp(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nphi,E_l(:,:,:,taup1,m)) 
   call setcyclic_xyp(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nphi,E_l(:,:,:,taup1,m))

   tauloc = taup1
  enddo

  !---------------------------------------------------------------------------------
  ! Diagnose Delta E
  !---------------------------------------------------------------------------------
  G_l(:,:,:,taup1,m) = E_l(:,:,:,taup1,m)*lambda_l(:,:,:,m)/(2*sign(abs(gamma_l(:,:,:,m))+1d-18,gamma_l(:,:,:,m)) )
  G_s(:,:,:,taup1,m) = E_s(:,:,:,taup1,m)*lambda_s(:,:,:,m)/(2*sign(abs(gamma_s(:,:,:,m))+1d-18,gamma_s(:,:,:,m)) )

  !if (rossmix_symmetrisation > 0d0 ) then
  ! G_l(:,:,:,taup1,m) = E_l(:,:,:,taup1,m)*2*gamma_l(:,:,:,m)/rossmix_symmetrisation
  ! G_s(:,:,:,taup1,m) = E_s(:,:,:,taup1,m)*2*gamma_s(:,:,:,m)/rossmix_symmetrisation
  !endif

  !if (rossmix_c_tau > 0d0 ) then
  ! do i=is_pe-onx,ie_pe+onx
  !  do j=js_pe-onx,je_pe+onx
  !   fxa = rossmix_c_tau*beta(i,j)*Rn(i,j,m)
  !   G_l(i,j,:,taup1,m) = E_l(i,j,:,taup1,m)*2*gamma_l(i,j,:,m)/fxa 
  !   G_s(i,j,:,taup1,m) = E_s(i,j,:,taup1,m)*2*gamma_s(i,j,:,m)/fxa
  !  enddo  
  ! enddo  
  !endif


 enddo


 !---------------------------------------------------------------------------------
 ! add effect of inverse energy cascade
 !---------------------------------------------------------------------------------

 if ( rossmix_inverse_energy > 0d0 ) then
   do m=2,nmodes
    ! E+ = E- - dt*a*E-*E+  ->  E+ = E-/(1+ dt*a*E-)
    aloc(:,:,:) = E_s(:,:,:,taup1,m)
    E_s(:,:,:,taup1,m) = E_s(:,:,:,taup1,m)/(1+dt_tracer*rossmix_inverse_energy*E_s(:,:,:,taup1,m))
    E_l(:,:,:,taup1,m) = E_l(:,:,:,taup1,m)+dt_tracer*rossmix_inverse_energy*E_s(:,:,:,taup1,m)*aloc(:,:,:)
   enddo
 endif

 !---------------------------------------------------------------------------------
 ! add effect of barotropisation
 !---------------------------------------------------------------------------------

 if ( rossmix_barotropisation > 0d0 ) then
   do m=nmodes,2,-1
    aloc(:,:,:) = E_s(:,:,:,taup1,m)
    E_s(:,:,:,taup1,m) = E_s(:,:,:,taup1,m)/(1+dt_tracer*rossmix_barotropisation*E_s(:,:,:,taup1,m))
    E_s(:,:,:,taup1,m-1) = E_s(:,:,:,taup1,m-1)+dt_tracer*rossmix_barotropisation*E_s(:,:,:,taup1,m)*aloc(:,:,:)

    aloc(:,:,:) = E_l(:,:,:,taup1,m)
    E_l(:,:,:,taup1,m) = E_l(:,:,:,taup1,m)/(1+dt_tracer*rossmix_barotropisation*E_l(:,:,:,taup1,m))
    E_s(:,:,:,taup1,m-1) = E_s(:,:,:,taup1,m-1)+dt_tracer*rossmix_barotropisation*E_l(:,:,:,taup1,m)*aloc(:,:,:)
   enddo
 endif


 
end subroutine integrate_rossmix







