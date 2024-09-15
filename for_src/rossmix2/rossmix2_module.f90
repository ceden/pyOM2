
module rossmix2_module   
      implicit none     
      logical :: enable_rossmix2 = .false.      
      logical :: enable_rossmix2_write_ub = .false.
      logical :: enable_rossmix2_simple   = .false.
      logical :: enable_rossmix2_lapack   = .false.
      logical :: enable_rossmix2_include_lateral_friction   = .false.
      logical :: enable_rossmix2_linear_damping  = .false.
      logical :: enable_rossmix2_appr_wavenumber = .false.
      logical :: enable_rossmix2_smooth_mom_fluxes = .true.
      logical :: enable_rossmix2_limit_energy = .false.
      logical :: enable_rossmix2_calc_phidot = .false.
      
      integer :: nphi = 10  ! grid points in angle direction
      real*8  :: rossmix2_fmin=0.05e-4
      real*8  :: rossmix2_N2min= (10*0.05e-4)**2
      real*8  :: rossmix2_calc_parameter_period = 86400. 
      real*8  :: rossmix2_lateral_friction_fraction = 0.5  ! was 0.25 in first attempt
      
      
      real*8  :: rossmix2_Kh = 1000.0   ! diffusivity for E in zonal/meridional direction
      real*8  :: rossmix2_Kp = 0.2e-6   ! same but in phi direction
      real*8  :: rossmix2_damping = 1./(30.*86400.) ! linear damping of E 
      real*8  :: rossmix2_ub_limiter = 5000.
      real*8  :: rossmix2_uu_limiter = 0.05
      real*8  :: rossmix2_mindepth = 800.0
      
      real*8 :: rossmix2_kh_appr_max = 4.0
      real*8 :: rossmix2_kh_appr_min = 0.0
      
      real*8 :: rossmix2_A1 = 0.340
      real*8 :: rossmix2_A2 = 0.794
      real*8 :: rossmix2_A3 = 0.138
      real*8 :: rossmix2_A4 = 0.206
      real*8 :: rossmix2_A7 = 0.2
      real*8 :: rossmix2_A8 = 0.0
      
      real*8 :: rossmix2_congr_epsilon = 1e-4
      
      real*8, allocatable, dimension(:,:,:)   :: maskTp     ! mask (1 or 0) for tracer points in (x,y,p) space
      real*8, allocatable, dimension(:,:,:)   :: maskUp     ! mask for U-grid points, where U is zonal group velocity
      real*8, allocatable, dimension(:,:,:)   :: maskVp     ! mask for V-grid points, where V is meridional group velocity
      real*8, allocatable, dimension(:,:,:)   :: maskWp     ! mask for W-grid points, where W is change of angle phi
      real*8, allocatable, dimension(:)       :: phit,dphit ! angle coordinate of T-grid 
      real*8, allocatable, dimension(:)       :: phiu,dphiu ! angle coordinate of U-grid 

      real*8, allocatable, dimension(:,:,:)   :: Nsqr_lim   ! limited Nsqr on W grid
      real*8, allocatable, dimension(:,:,:)   :: cgu,cgv,phidot  ! group velocity   
      !real*8, allocatable, dimension(:,:)     :: A1,A2      ! spectral shape parameter
      real*8, allocatable, dimension(:,:,:,:) :: ub,uu      ! eddy flux structure
      real*8, allocatable, dimension(:,:,:,:) :: E_r        ! wave energy 
      complex*16, allocatable, dimension(:,:,:)   :: om_max ! frequency 
      real*8, allocatable, dimension(:,:,:)   :: Bx,By      ! eddy streamfunction
      real*8, allocatable, dimension(:,:,:)   :: ue,ve,we   ! eddy velocities
      real*8, allocatable, dimension(:,:,:)   :: ug,vg      ! background velocities
      real*8, allocatable, dimension(:,:,:)   :: forcx,forcy! energy transfers
      real*8, allocatable, dimension(:,:,:)   :: forcuu
      real*8, allocatable, dimension(:,:)     :: forc_diss
      real*8, allocatable, dimension(:,:,:)   :: stress_ux,stress_uy
      real*8, allocatable, dimension(:,:,:)   :: stress_vx,stress_vy
      real*8, allocatable, dimension(:,:)     :: U1,V1
      real*8, allocatable, dimension(:,:)     :: kh_appr

      real*8, allocatable, dimension(:,:,:)   :: ub_lim, uu_lim

      real*8, allocatable, dimension(:,:)     :: R1         ! 1. Rossby radius 
      real*8, allocatable, dimension(:,:,:)   :: phi1       ! 1. vertical mode 
      real*8, allocatable, dimension(:,:,:)   :: phi1W      ! 1. vertical mode, but on W grid 
      real*8, allocatable, dimension(:,:,:)   :: phiz1      ! derivative of 1. vertical mode 
      real*8, allocatable, dimension(:,:,:)   :: phiz1W     ! derivative of 1. vertical mode, but on W grid 
      
end module rossmix2_module   



subroutine allocate_rossmix2_module
 use main_module   
 use rossmix2_module   
 implicit none
 
 if (.not. enable_rossmix2) return
 
 allocate( phit(nphi), dphit(nphi), phiu(nphi), dphiu(nphi))  
 phit=0;dphit=0;phiu=0;dphiu=0
 allocate( maskTp(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nphi) )
 allocate( maskUp(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nphi) )
 allocate( maskVp(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nphi) )
 allocate( maskWp(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nphi) )
 maskTp=0; maskUp=0; maskVp=0; maskWp=0; 

 allocate( Nsqr_lim(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz)); Nsqr_lim = 0.0
 allocate( cgu(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nphi) ); cgu=0
 allocate( cgv(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nphi) ); cgv=0
 allocate( phidot(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nphi) ); phidot=0
 !allocate( A1(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx) ); A1=0
 !allocate( A2(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx) ); A2=0
 allocate( ub(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz,nphi)); ub = 0.0
 allocate( uu(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz,nphi)); uu = 0.0
 allocate( E_r(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nphi,3)); E_r= 0.0 
 allocate( om_max(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nphi)); om_max= 0.0 
 allocate( Bx(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz)); Bx = 0.0
 allocate( By(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz)); By = 0.0 
 allocate( ue(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz)); ue = 0.0
 allocate( ve(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz)); ve = 0.0
 allocate( we(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz)); we = 0.0
 allocate( ug(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz)); ug = 0.0
 allocate( vg(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz)); vg = 0.0
 allocate( forcx(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nphi) ); forcx=0
 allocate( forcy(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nphi) ); forcy=0
 allocate( forcuu(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nphi) ); forcuu=0
 allocate( forc_diss(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx) ); forc_diss=0

 allocate( stress_ux(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); stress_ux=0
 allocate( stress_uy(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); stress_uy=0
 allocate( stress_vx(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); stress_vx=0
 allocate( stress_vy(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); stress_vy=0

 allocate(U1(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx)); U1=0
 allocate(V1(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx)); V1=0
 allocate(kh_appr(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx)); kh_appr=0

 allocate( ub_lim(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); ub_lim=0
 allocate( uu_lim(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); uu_lim=0
 
 allocate( R1(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx) ); R1 = 0.0
 allocate( phi1(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); phi1 = 0.0
 allocate( phi1W(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); phi1W = 0.0
 allocate( phiz1(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); phiz1 = 0.0
 allocate( phiz1W(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); phiz1W = 0.0
end subroutine allocate_rossmix2_module