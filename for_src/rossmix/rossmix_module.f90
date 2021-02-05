

module rossmix_module   
      ! main module containing all relevant arrays and parameter
      implicit none

      integer :: nphi     ! grid points in angle direction
      integer :: nmodes   ! number of vertical modes


      logical :: enable_rossmix = .false.
      logical :: enable_rossmix_mean_flow_interaction = .true.
      logical :: enable_rossmix_vertical_stress       = .false.
      logical :: enable_rossmix_lateral_stress        = .true.
      logical :: enable_rossmix_diffusion             = .false.
      logical :: enable_rossmix_use_gamma             = .false.


      real*8  :: rossmix_calc_cg_int = 86400.0              ! interval in sec. to calculate group velocties, etc.
      integer :: rossmix_barotropic_fac = 1                 ! reduction factor for barotropic time stepping
      real*8  :: rossmix_inverse_energy  = 0.0              ! parameter for inverse energy cascade
      real*8  :: rossmix_barotropisation = 0.0              ! parameter for inverse energy cascade
      real*8  :: rossmix_symmetrisation  = 0.0              ! inverse time scale for vertical symmetrisation
      real*8  :: rossmix_c_tau    = 0.0 
      real*8  :: rossmix_c_lambda = 1.0 

      real*8  :: rossmix_rayleigh_friction = 0.0            ! rayleigh friction coefficient
      real*8  :: rossmix_lat_diffusion = 0.0                ! lateral diffusivity for wave energy
      real*8  :: rossmix_phi_diffusion = 0.0                ! lateral diffusivity for wave energy
      real*8  :: rossmix_fmin  = 1e-5                       ! minimal value of Coriolis parameter
      real*8  :: rossmix_N2min = 1e-7                       ! minimal value of Coriolis parameter
      real*8  :: rossmix_hmin = 200.0

      real*8, allocatable, dimension(:,:,:)   :: maskTp     ! mask (1 or 0) for tracer points in (x,y,p) space
      real*8, allocatable, dimension(:,:,:)   :: maskUp     ! mask for U-grid points, where U is zonal group velocity
      real*8, allocatable, dimension(:,:,:)   :: maskVp     ! mask for V-grid points, where V is meridional group velocity
      real*8, allocatable, dimension(:,:,:)   :: maskWp     ! mask for W-grid points, where W is change of angle phi

      real*8, allocatable, dimension(:)       :: phit,dphit ! angle coordinate of T-grid 
      real*8, allocatable, dimension(:)       :: phiu,dphiu ! angle coordinate of U-grid 

      real*8, allocatable, dimension(:,:,:)   :: Rn          ! Rossby radius for each mode in m
      real*8, allocatable, dimension(:,:,:,:) :: phin,phizn  ! eigenfunction on T grid points
      real*8, allocatable, dimension(:,:,:,:) :: phinW,phiznW! eigenfunction on W grid points
      real*8, allocatable, dimension(:,:,:,:) :: phi_ga_l    ! structure function
      real*8, allocatable, dimension(:,:,:,:) :: phi_ga_s    ! structure function
      real*8, allocatable, dimension(:,:,:)   :: phi_nu_l    ! structure function
      real*8, allocatable, dimension(:,:,:)   :: phi_nu_s    ! structure function

      real*8, allocatable, dimension(:,:,:) ::   ray_fric_s ! for rayleigh friction
      real*8, allocatable, dimension(:,:,:) ::   ray_fric_l ! for rayleigh friction

      real*8, allocatable, dimension(:,:,:,:,:) :: E_s,E_l   ! sum wave energy per angle in J/m^2 (or m^3/s^2 )
      real*8, allocatable, dimension(:,:,:,:,:) :: G_s,G_l   ! diff in wave energy per angle in J/m^2 (or m^3/s^2 )

      real*8, allocatable                     :: cgu_s(:,:,:,:)        
      real*8, allocatable                     :: cgv_s(:,:,:,:)        
      real*8, allocatable                     :: phidot_s(:,:,:,:)        
      real*8, allocatable                     :: cgu_l(:,:,:,:)        
      real*8, allocatable                     :: cgv_l(:,:,:,:)        
      real*8, allocatable                     :: phidot_l(:,:,:,:)        

      real*8, allocatable                     :: gamma_l(:,:,:,:)        
      real*8, allocatable                     :: gamma_s(:,:,:,:)        
      real*8, allocatable                     :: lambda_l(:,:,:,:)        
      real*8, allocatable                     :: lambda_s(:,:,:,:)        
      real*8, allocatable                     :: nu_s(:,:,:,:)        
      real*8, allocatable                     :: nu_l(:,:,:,:)        

      real*8, allocatable                     :: Bx(:,:,:)
      real*8, allocatable                     :: By(:,:,:)
      real*8, allocatable                     :: stress_ux(:,:,:)        
      real*8, allocatable                     :: stress_uy(:,:,:)        
      real*8, allocatable                     :: stress_vx(:,:,:)        
      real*8, allocatable                     :: stress_vy(:,:,:)        

      real*8, allocatable                     :: cfac1_s(:,:,:),cfac1_l(:,:,:)
      real*8, allocatable                     :: cfac2_s(:,:,:),cfac2_l(:,:,:)
      real*8, allocatable                     :: km_l(:,:,:)
      real*8, allocatable                     :: rossmix_ue(:,:,:),rossmix_ve(:,:,:),rossmix_we(:,:,:)
      real*8, allocatable                     :: P_diss_ross(:,:,:)
      real*8, allocatable                     :: K_diss_ross(:,:,:)

      integer, allocatable                     :: boundary_west_i(:),boundary_west_j(:)
      integer, allocatable                     :: boundary_east_i(:),boundary_east_j(:)
      integer, allocatable                     :: boundary_north_i(:),boundary_north_j(:)
      integer, allocatable                     :: boundary_south_i(:),boundary_south_j(:)

end module rossmix_module   


subroutine allocate_rossmix_module
 ! allocate all arrays within main module
 use main_module   
 use rossmix_module   
 implicit none

 if (enable_rossmix) then

  allocate( phit(nphi), dphit(nphi), phiu(nphi), dphiu(nphi)) ; phit=0;dphit=0;phiu=0;dphiu=0

  allocate( maskTp(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nphi) )
  allocate( maskUp(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nphi) )
  allocate( maskVp(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nphi) )
  allocate( maskWp(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nphi) )
  maskTp=0; maskUp=0; maskVp=0; maskWp=0; 

  allocate( Rn(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nmodes) ); Rn = 0.0
  allocate( phin(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz,nmodes) ); phin = 0.0
  allocate( phinW(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz,nmodes) ); phinW = 0.0
  allocate( phizn(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz,nmodes) ); phizn = 0.0
  allocate( phiznW(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz,nmodes) ); phiznW = 0.0


  allocate( phi_ga_l(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz,nmodes) ); phi_ga_l = 0.0
  allocate( phi_ga_s(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz,nmodes) ); phi_ga_s = 0.0
  allocate( phi_nu_l(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nmodes) ); phi_nu_l = 0.0
  allocate( phi_nu_s(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nmodes) ); phi_nu_s = 0.0


  allocate( E_s(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nphi,3,nmodes) ); E_s=0
  allocate( E_l(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nphi,3,nmodes) ); E_l=0

  allocate( G_s(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nphi,3,nmodes) ); G_s=0
  allocate( G_l(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nphi,3,nmodes) ); G_l=0

  allocate( cgu_s(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nphi,nmodes) ); cgu_s=0
  allocate( cgv_s(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nphi,nmodes) ); cgv_s=0
  allocate( phidot_s(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nphi,nmodes) ); phidot_s=0

  allocate( cgu_l(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nphi,nmodes) ); cgu_l=0
  allocate( cgv_l(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nphi,nmodes) ); cgv_l=0
  allocate( phidot_l(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nphi,nmodes) ); phidot_l=0

  allocate( gamma_l(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nphi,nmodes) ); gamma_l=0
  allocate( gamma_s(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nphi,nmodes) ); gamma_s=0
  allocate( nu_l(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nphi,nmodes) ); nu_l=0
  allocate( nu_s(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nphi,nmodes) ); nu_s=0
  allocate( lambda_l(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nphi,nmodes) ); lambda_l=0
  allocate( lambda_s(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nphi,nmodes) ); lambda_s=0

  allocate( Bx(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); Bx=0
  allocate( By(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); By=0

  allocate( stress_ux(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); stress_ux=0
  allocate( stress_uy(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); stress_uy=0
  allocate( stress_vx(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); stress_vx=0
  allocate( stress_vy(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); stress_vy=0

  allocate( cfac1_s(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nmodes), cfac1_l(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nmodes) )
  cfac1_s = 0; cfac1_l = 0
  allocate( cfac2_s(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nmodes), cfac2_l(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nmodes) )
  cfac2_s = 0; cfac2_l = 0
  allocate( km_l(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nmodes) ); km_l = 0

  allocate( ray_fric_s(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nmodes))
  allocate( ray_fric_l(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nmodes))
  ray_fric_s = 0; ray_fric_l = 0

  allocate( rossmix_ue(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); rossmix_ue=0
  allocate( rossmix_ve(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); rossmix_ve=0
  allocate( rossmix_we(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); rossmix_we=0

  allocate( P_diss_ross(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); P_diss_ross=0
  allocate( K_diss_ross(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); K_diss_ross=0
 endif
end subroutine allocate_rossmix_module
