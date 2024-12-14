

module biharmonic_thickness_module   
!=======================================================================
!=======================================================================
  implicit none
  logical :: biharmonic_thickness_active = .false.
  logical :: enable_biharmonic_thickness_explicit = .false. ! enable biharmonic thickness mixing (in res. mom. formulation)
  logical :: enable_biharmonic_thickness_bolus_skew = .false. ! enable biharmonic thickness mixing, bolus formulation skew  
  logical :: enable_biharmonic_thickness_backscatter_harmonic    = .false.
  logical :: enable_biharmonic_thickness_backscatter_skew = .false.
  logical :: enable_biharmonic_thickness_backscatter_integrate_energy  = .false.
  logical :: enable_biharmonic_thickness_backscatter_taper_mld = .false.
  integer :: biharmonic_thickness_backscatter_smooth    = 1
  
  real*8  :: biha_dslope=1d-9
  real*8  :: biha_slopec=2d-9 
  real*8  :: A_thk_0 = 0.
  real*8  :: mld0_thk = 75.0
  real*8  :: e_back_diff = 1000.0
  real*8  :: A_thk_max = 500.
  real*8  :: L_back = 2e3
  
  logical :: enable_biharmonic_thickness_cut_upper = .false.         ! cut upper XX meter for explicit formulation
  logical :: enable_biharmonic_thickness_scale_A_thkbi_cut = .false. ! scale that parameter with cos(lat)**2 for explicit 

  real*8 :: A_thkbi = 0.0                   ! biharmonic thickness diffusivity in m^4/s
  real*8 :: thkbi_f_min = 1e-6              ! threshold for Coriolis parameter in expl. bihar. thickness formulation in 1/s
  real*8 :: K_thkbi_clip = 5.               ! clipping for visc. in multiples of f, only biharmonic_thickness_mixing_explicit
      
  integer:: biharmonic_thickness_mixing_iter = 1 ! iterations in sub loop for the explicit formulation
  real*8 :: A_thkbi_cut = 0.5
  real*8 :: A_thkbi_cut_upper_len = 35.

  real*8,allocatable :: K_thkbi(:,:,:)    ! actual A_thkbi after clipping
  real*8,allocatable :: A_thk(:,:,:)      ! 
  real*8,allocatable :: diss_biha_skew(:,:,:)  
  real*8,allocatable :: diss_back(:,:,:) 
  
  real*8,allocatable :: mld(:,:),mldu(:,:),mldv(:,:)
  real*8,allocatable :: e_back(:,:,:)
  
  real*8,allocatable :: diag_Bx(:,:,:),diag_By(:,:,:)
  real*8,allocatable :: diag_Bx_back(:,:,:),diag_By_back(:,:,:)
end module biharmonic_thickness_module  





subroutine allocate_biharmonic_thickness_module
  use main_module
  use biharmonic_thickness_module 
  implicit none
      
  if ( .not. biharmonic_thickness_active) return    
        
  allocate( K_thkbi(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); K_thkbi = 0.0
  allocate( mld(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx) ) ; mld=0d0
  allocate( mldu(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx) ); mldu=0d0
  allocate( mldv(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx) ); mldv=0d0
  
  allocate( diag_Bx(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); diag_Bx=0d0
  allocate( diag_By(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); diag_By=0d0

  if (   enable_biharmonic_thickness_backscatter_skew &
    .or. enable_biharmonic_thickness_backscatter_harmonic ) then
   allocate( A_thk(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); A_thk = 0.0
   allocate( diag_Bx_back(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); diag_Bx_back = 0.0
   allocate( diag_By_back(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); diag_By_back = 0.0
   if (enable_conserve_energy) then 
     allocate( diss_back(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); diss_back = 0
   endif     
   if (enable_biharmonic_thickness_backscatter_integrate_energy) then
     allocate( e_back(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,3) ); e_back=0d0
   endif
  endif
  
  if ( enable_biharmonic_thickness_bolus_skew  .and. enable_conserve_energy) then
   allocate( diss_biha_skew(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); diss_biha_skew = 0
  endif
  

end subroutine allocate_biharmonic_thickness_module

