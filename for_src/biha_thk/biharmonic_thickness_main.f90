

subroutine biharmonic_thickness_friction
 use biharmonic_thickness_module
 implicit none  
 
 if (enable_biharmonic_thickness_bolus_skew .or. enable_biharmonic_thickness_backscatter_skew &
     .or. enable_biharmonic_thickness_backscatter_harmonic) then
   call biha_find_ml_depth
 endif
 
 if (enable_biharmonic_thickness_explicit) call biharmonic_thickness_explicit
 if (enable_biharmonic_thickness_backscatter_harmonic) call biharmonic_thickness_backscatter_harmonic
end subroutine biharmonic_thickness_friction


subroutine biharmonic_thickness_mixing
 use main_module 
 use biharmonic_thickness_module
 implicit none  
 if (enable_biharmonic_thickness_bolus_skew ) call biharmonic_thickness_skew_tensor
end subroutine biharmonic_thickness_mixing




subroutine biharmonic_thickness_initialize
 use main_module
 use biharmonic_thickness_module 
 implicit none
 real*8 :: fxa,dt_loc,fxb
 integer :: i,j,k

 if   ( enable_biharmonic_thickness_explicit &
   .or. enable_biharmonic_thickness_bolus_skew &
   .or. enable_biharmonic_thickness_backscatter_skew &
   .or. enable_biharmonic_thickness_backscatter_harmonic) then
    biharmonic_thickness_active = .true.
 else
    biharmonic_thickness_active = .false.
 endif
 
 if ( .not. biharmonic_thickness_active) return

 if (enable_biharmonic_thickness_backscatter_skew .and. &
     .not. enable_biharmonic_thickness_bolus_skew) then
   if (my_pe==0) print*,'ERROR (1) in biharmonic_thickness_initialize, check source'
   call halt_stop('in biharmonic_thickness_initialize (1) ')
 endif    

 if (enable_biharmonic_thickness_backscatter_integrate_energy .and. &
     .not. enable_conserve_energy) then
   if (my_pe==0) print*,'ERROR (2) in biharmonic_thickness_initialize, check source'
   call halt_stop('in biharmonic_thickness_initialize (2) ')  
 endif

 
 if (my_pe==0) print'(a)',' setting up biharmonic thickness mixing'


 call allocate_biharmonic_thickness_module
  
 call init_snap_biharmonic_thickness

 if (enable_biharmonic_thickness_explicit) then
   fxa = 0.
   dt_loc = dt_mom/biharmonic_thickness_mixing_iter
   do k=1,nz
    do j=js_pe,je_pe
     do i=is_pe,ie_pe
      fxb = A_thkbi*(cost(j)**biha_friction_cosPower)**2 * (1./10.)**2
      fxa = max( fxa, fxb*dt_loc*2/(cost(j)*dxt(i))**2*2/dzt(k)**2 )
      fxa = max( fxa, fxb*dt_loc*2/(        dyt(j))**2*2/dzt(k)**2 )
     enddo 
    enddo
   enddo 
   call global_max(fxa) 
   if (my_pe==0) print*,' biharm. thickness diffusion max criterium = ',fxa
   if (fxa>1) then
     if (my_pe==0) print*,'ERROR: A_thkbi too large,  criterium must be <1'
     !call halt_stop('in setup')
   endif
 endif
  

end subroutine biharmonic_thickness_initialize




subroutine biha_find_ml_depth!(mld,mldu,mldv)
 use main_module    
 use biharmonic_thickness_module
 implicit none
 integer :: i,j,k
 real*8 :: N2_min,N2_int,fxa,get_rho
 integer :: ml_algorithm =1
 real*8  :: ml_Nsqr_fac   = 8d0  ! for mixed layer criterion
 real*8  :: ml_Delta_rho  = 0.02 ! for mixed layer criterion

  !-----------------------------------------------------------------------
  ! find mixed layer depth 
  !-----------------------------------------------------------------------
  if (ml_algorithm == 1) then
    ! as given in Fox-Kemper et al. but with modification wrt isoneutral_ml
    mld=dzw(nz)/2.
    do j=js_pe-onx,je_pe+onx
     do i=is_pe-onx,ie_pe+onx      
       if (kbot(i,j)>0 .and. Nsqr(i,j,nz,tau) < 1d-12) then   
         N2_min = max(0d0,Nsqr(i,j,nz,tau)); N2_int=0d0; 
         do k=nz-1,kbot(i,j),-1
            fxa = max(0d0,Nsqr(i,j,k,tau))
            if ( fxa-N2_min  .gt. ml_Nsqr_fac*N2_int/mld(i,j)) exit
            N2_min = min(fxa,N2_min)
            N2_int = N2_int + fxa*dzw(k)
            mld(i,j)  = mld(i,j) + dzw(k)
         enddo
       endif
     enddo
    enddo
   
 elseif (ml_algorithm == 3) then
    do j=js_pe-onx,je_pe+onx
     do i=is_pe-onx,ie_pe+onx
         if (kbot(i,j)>0) then
          do k=nz-1,kbot(i,j),-1
            mld(i,j) = zw(nz)-zw(k)
            if ( get_rho(salt(i,j,k,tau),temp(i,j,k,tau),abs(zt(nz)))  .gt.  rho(i,j,nz,tau) +ml_Delta_rho ) exit
          enddo
         endif
     enddo
    enddo
  else
    call halt_stop('ERROR in biha_find_ml_depth: wrong parameter for ml_algorithm ')
  endif

  mldu=mld; mldv=mld
  do i=is_pe-onx,ie_pe+onx-1
    mldu(i,:) = 0.5*(mld(i+1,:)+mld(i,:)) 
  enddo
  do j=js_pe-onx,je_pe+onx-1
    mldv(:,j) = 0.5*(mld(:,j+1)+mld(:,j)) 
  enddo
end subroutine biha_find_ml_depth




