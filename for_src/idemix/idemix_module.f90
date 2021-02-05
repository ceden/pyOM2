
module idemix_module   
!=======================================================================
! module containing all relevant arrays and parameter for IDEMIX
!=======================================================================
      implicit none

!---------------------------------------------------------------------------------
!     Idemix 1.0
!---------------------------------------------------------------------------------
      logical :: enable_idemix = .false.
      real*8, allocatable :: dE_iw(:,:,:,:)                   ! tendency due to advection using Adam Bashforth
      real*8, allocatable :: E_iw(:,:,:,:)                    ! total energy of internal waves
      real*8, allocatable :: c0(:,:,:)                        ! mean vertical group velocity
      real*8, allocatable :: v0(:,:,:)                        ! mean lateral group velocity
      real*8, allocatable :: cstar(:,:)                       ! modal gravity wave speed
      real*8, allocatable :: alpha_c(:,:,:)                   ! dissipation parameter
      real*8, allocatable :: forc_iw_bottom(:,:)              ! energy flux into IW at bottom
      real*8, allocatable :: forc_iw_surface(:,:)             ! energy flux into IW at surface
      real*8, allocatable :: iw_diss(:,:,:)                   ! total dissipation of IW energy
      real*8 :: tau_v=1.0*86400.0                             ! time scale for vertical symmetrisation
      real*8 :: tau_h=15.0*86400.0                            ! time scale for horizontal symmetrisation
      real*8 :: gamma=1.57                                    ! 
      real*8 :: jstar = 10.0                                  ! spectral bandwidth in modes
      real*8 :: mu0   = 4.0/3.0                               ! dissipation parameter
      logical :: enable_idemix_hor_diffusion = .false.        ! lateral diffusion in IDEMIX 1.0
      logical :: enable_idemix_hor_diffusion_iter = .false.   ! lateral diffusion with iterative method in case of high resolution
      integer :: idemix_hor_diffusion_iter = 1                ! number of iterations for this case
      logical :: enable_eke_diss_bottom = .false.             ! dissipate of EKE is injected at the bottom into internal wave field
      logical :: enable_eke_diss_surfbot = .false.            ! same but EKE is injected at surface and bottom
      real*8  :: eke_diss_surfbot_frac = 1.0                  ! fraction which goes into bottom
      logical :: enable_idemix_superbee_advection   = .false. ! use superbee advection for internal wave energy
      logical :: enable_idemix_dst3_advection       = .false. ! use 2. order direct space-time advection for internal wave energy
      logical :: enable_idemix_AB_time_stepping     = .true.  ! use Adams-Bashforth time stepping for advection

      
!---------------------------------------------------------------------------------
!     IDEMIX 2.0
!---------------------------------------------------------------------------------
      logical :: enable_idemix_M2     = .false.
      logical :: enable_idemix_niw    = .false.

      integer :: np = 0
      real*8, allocatable :: phit(:),phiu(:),dphit(:),dphiu(:)
      real*8, allocatable :: maskTp(:,:,:)
      real*8, allocatable :: maskUp(:,:,:)
      real*8, allocatable :: maskVp(:,:,:)
      real*8, allocatable :: maskWp(:,:,:)
      integer, allocatable, dimension(:,:,:)  :: bc_south   ! index for southern reflection boundary condition
      integer, allocatable, dimension(:,:,:)  :: bc_north   ! same for northern boundaries
      integer, allocatable, dimension(:,:,:)  :: bc_east    ! same for eastern boundaries
      integer, allocatable, dimension(:,:,:)  :: bc_west    ! same for western boundaries

      real*8, allocatable :: topo_hrms(:,:),topo_lam(:,:),topo_shelf(:,:)
      real*8, allocatable :: E_M2(:,:,:,:),dE_M2p(:,:,:,:)
      real*8, allocatable :: E_niw(:,:,:,:),dE_niwp(:,:,:,:)

      real*8, allocatable :: tau_niw(:,:),forc_niw(:,:,:)
      real*8, allocatable :: tau_M2(:,:),forc_M2(:,:,:)
      real*8, allocatable :: alpha_M2_cont(:,:)
      real*8, allocatable :: M2_psi_diss(:,:,:)  ! dissipation by PSI 

      real*8, allocatable :: cn(:,:),phin(:,:,:),phinz(:,:,:)
      real*8, allocatable :: omega_niw(:,:)
      real*8   :: omega_M2 !=  2*pi/( 12*60*60 +  25.2 *60 )   ! M2 frequency in 1/s
      real*8, allocatable :: cg_niw(:,:),cg_M2(:,:)
      real*8, allocatable :: kdot_y_M2(:,:),kdot_y_niw(:,:)
      real*8, allocatable :: kdot_x_M2(:,:),kdot_x_niw(:,:)
      real*8, allocatable :: u_M2(:,:,:),v_M2(:,:,:),w_M2(:,:,:)
      real*8, allocatable :: u_niw(:,:,:),v_niw(:,:,:),w_niw(:,:,:)
      real*8, allocatable :: E_struct_niw(:,:,:),E_struct_M2(:,:,:)
      real*8, allocatable :: E_niw_int(:,:),E_M2_int(:,:)
!---------------------------------------------------------------------------------
!     IDEMIX 3.0
!---------------------------------------------------------------------------------
      real*8 :: mu0_min = 1d-4, tc_max = 1d-3, Noverf_min0=10., Noverf_min1=0.75
      logical :: enable_idemix3               = .false.
      real*8, allocatable :: viscU_idemix(:,:,:)     !  
      real*8, allocatable :: viscV_idemix(:,:,:)     !  
      real*8, allocatable :: K_diss_idemix(:,:,:)    !  
      real*8, allocatable :: c0_u(:,:,:),c0_v(:,:,:) !  
      real*8, allocatable :: v0_u(:,:,:),v0_v(:,:,:) !  
      real*8, allocatable :: alpha_c_u(:,:,:),alpha_c_v(:,:,:)   !  
      real*8, allocatable :: Tc_U(:,:,:),Tc_V(:,:,:) !  
      real*8, allocatable :: Om_id3_u(:,:,:),Om_id3_v(:,:,:)       !  
      real*8, allocatable :: G_d(:,:,:,:),G_s(:,:,:,:)   !  
      real*8, allocatable :: F_d(:,:,:,:),F_s(:,:,:,:)   !  
      real*8, allocatable :: O_d(:,:,:,:),O_s(:,:,:,:)   !  
      real*8, allocatable :: P_d(:,:,:,:),P_s(:,:,:,:)   !  
      real*8, allocatable :: forc_iw_bottom_u(:,:)              ! energy flux into IW at bottom
      real*8, allocatable :: forc_iw_surface_u(:,:)             ! energy flux into IW at surface
      real*8, allocatable :: forc_iw_bottom_v(:,:)              ! energy flux into IW at bottom
      real*8, allocatable :: forc_iw_surface_v(:,:)             ! energy flux into IW at surface
      logical :: enable_idemix3_obc_taper_south  = .false.   ! taper wave drag to zero near open boundary conditions
      logical :: enable_idemix3_obc_taper_north  = .false.
      real*8 :: id3_north_taper_len=2, id3_south_taper_len=2
      logical :: enable_idemix3_cut_negative_energy = .false.

!---------------------------------------------------------------------------------
!     Lee waves
!---------------------------------------------------------------------------------
      logical :: enable_leewaves          = .false. ! enables leewaves   
      logical :: enable_leewaves_obc_taper_north    ! taper leewaves to zero near open boundary 
      logical :: enable_leewaves_obc_taper_south 
      logical :: enable_leewaves_slowflux = .false.
      real*8 :: leewaves_north_taper_len=2
      real*8 :: leewaves_south_taper_len=2
      real*8 :: tau_v_lee = 3*86400.
      integer :: idemix_lee_ts_fac = 1              ! sub loop time stepping for lee waves
      integer :: nph = 32                           ! grid points in angle space.
      integer :: nk  = 32                           ! grid points in wavenumber space. 
      real*8              :: nu = 0.9               ! Hurst number
      real*8              :: inv_fr_c = 0.75        ! critical inverse Froude number
      real*8  ::  Ah_lee = 50.
      real*8  ::  u0_lee_min = 1d-4
      real*8, allocatable, dimension(:)    :: ph    !angle of propagation 
      real*8, allocatable, dimension(:,:)  :: lam   !Topographic wavelength
      real*8, allocatable, dimension(:,:)  :: k_s   !Topographic wavenumber
      real*8, allocatable, dimension(:,:)  :: h_rms !RMS-height of topographic spectrum
      real*8, allocatable, dimension(:,:)  :: k_n   !Topographic wavenumber (normal direction)
      real*8, allocatable, dimension(:,:)  :: ph_s  !strike angle (used in 12.338)
      
      real*8, allocatable :: flux_lee_tot(:,:)
      real*8, allocatable :: inv_fr(:,:)            ! inverse Froude number 
      real*8, allocatable :: u_bot(:,:)             ! Bottom velocity in u-direction
      real*8, allocatable :: v_bot(:,:)             ! bottom velocity in v-direction
      real*8, allocatable :: u0_bot(:,:)
      real*8, allocatable :: topo_fac(:,:)
      real*8, allocatable :: c_lee(:,:,:),Tc_lee(:,:,:) ,Lambda_0(:,:,:)
      real*8, allocatable :: E_lee_d(:,:,:,:),dE_lee_d(:,:,:,:) 
      real*8, allocatable :: E_lee_s(:,:,:,:),dE_lee_s(:,:,:,:) 
      real*8, allocatable :: E_lee_p(:,:,:,:),dE_lee_p(:,:,:,:)
      real*8, allocatable :: E_lee_m(:,:,:,:),dE_lee_m(:,:,:,:) 
      real*8, allocatable :: mean_flow_to_E_lee(:,:,:) 
      real*8, allocatable :: iw_diss_lee(:,:,:) 
      real*8, allocatable :: stress_lee_u(:,:),stress_lee_v(:,:)
      
end module idemix_module   



subroutine allocate_idemix_module
  use main_module
  use idemix_module
  implicit none

  if (enable_idemix) then
   allocate(dE_iw(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz,3) );dE_iw = 0
   allocate( E_iw(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz,3) ); E_iw = 0
   allocate( cstar(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx) ); cstar = 0
   allocate( c0(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); c0 = 0
   allocate( v0(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); v0 = 0
   allocate( alpha_c(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); alpha_c = 0
   allocate( iw_diss(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); iw_diss = 0
   allocate( forc_iw_surface(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx) ); forc_iw_surface = 0
   allocate( forc_iw_bottom(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx) ); forc_iw_bottom = 0
   
   if (enable_idemix3) then
     allocate( viscU_idemix(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); viscU_idemix = 0
     allocate( viscV_idemix(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); viscV_idemix = 0
     allocate( K_diss_idemix(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); K_diss_idemix = 0
     allocate( c0_u(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); c0_u = 0
     allocate( c0_v(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); c0_v = 0
     allocate( v0_u(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); v0_u = 0
     allocate( v0_v(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); v0_v = 0
     allocate( alpha_c_u(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); alpha_c_u = 0
     allocate( alpha_c_v(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); alpha_c_v = 0

     allocate( Om_id3_u(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); Om_id3_u = 0
     allocate( Om_id3_v(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); Om_id3_v = 0
     allocate( Tc_U(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); Tc_U = 0
     allocate( Tc_V(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); Tc_V = 0
     
     allocate( G_d(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz,3) ); G_d = 0
     allocate( F_d(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz,3) ); F_d = 0
     allocate( G_s(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz,3) ); G_s = 0
     allocate( F_s(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz,3) ); F_s = 0
     allocate( O_d(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz,3) ); O_d = 0
     allocate( P_d(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz,3) ); P_d = 0
     allocate( O_s(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz,3) ); O_s = 0
     allocate( P_s(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz,3) ); P_s = 0
     allocate( forc_iw_surface_u(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx) ); forc_iw_surface_u = 0
     allocate( forc_iw_bottom_u(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx) ); forc_iw_bottom_u = 0
     allocate( forc_iw_surface_v(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx) ); forc_iw_surface_v = 0
     allocate( forc_iw_bottom_v(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx) ); forc_iw_bottom_v = 0
   endif


   if (enable_leewaves) then

     allocate( ph(nph) );   ph = 0.                                     
     allocate( k_s(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx) );       k_s   = 0.0
     allocate( h_rms(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx) );     h_rms = 0.0
     allocate( k_n(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx) );       k_n   = 0.0
     allocate( ph_s(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx) );      ph_s  = 0.0
     allocate( flux_lee_tot(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx) ); flux_lee_tot = 0
     allocate( inv_fr(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx) ); inv_fr = 0
     allocate( u_bot(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx) ); u_bot = 0
     allocate( v_bot(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx) ); v_bot = 0
     allocate( u0_bot(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx) ); u0_bot = 0
     allocate( stress_lee_u(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx) ); stress_lee_u = 0
     allocate( stress_lee_v(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx) ); stress_lee_v = 0
     allocate( topo_fac(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx) ); topo_fac = 0
     
     allocate(  E_lee_d(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz,3) );  E_lee_d = 0
     allocate( dE_lee_d(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz,3) ); dE_lee_d = 0
     allocate(  E_lee_s(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz,3) );  E_lee_s = 0
     allocate( dE_lee_s(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz,3) ); dE_lee_s = 0
     
     allocate(  E_lee_p(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz,3) );  E_lee_p = 0
     allocate(  E_lee_m(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz,3) );  E_lee_m = 0
     allocate( dE_lee_p(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz,3) ); dE_lee_p = 0
     allocate( dE_lee_m(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz,3) ); dE_lee_m = 0

     allocate( Lambda_0(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); Lambda_0 = 0
     allocate( c_lee(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); c_lee = 0
     allocate( Tc_lee(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); Tc_lee = 0
     allocate( mean_flow_to_E_lee(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); mean_flow_to_E_lee = 0
     allocate( iw_diss_lee(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); iw_diss_lee = 0
   endif
  endif


  if (enable_idemix_M2 .or. enable_idemix_niw ) then

    allocate( topo_shelf(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx) ); topo_shelf = 0
    allocate( topo_hrms(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx) ); topo_hrms = 0
    allocate( topo_lam(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx) ); topo_lam = 0
    allocate( phit(np), dphit(np), phiu(np), dphiu(np)) ; phit=0;dphit=0;phiu=0;dphiu=0
    allocate( maskTp(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,np) ); maskTp = 0
    allocate( maskUp(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,np) ); maskUp = 0
    allocate( maskVp(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,np) ); maskVp = 0
    allocate( maskWp(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,np) ); maskWp = 0
    allocate( cn(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx) ); cn = 0.0
    allocate( phin(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); phin = 0.0
    allocate( phinz(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); phinz = 0.0
    allocate( tau_M2(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx) ); tau_M2 = 0
    allocate( tau_niw(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx) ); tau_niw = 0
    allocate( alpha_M2_cont(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx) ); alpha_M2_cont = 0
    allocate( bc_south(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,np) )
    allocate( bc_north(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,np) )
    allocate( bc_west (is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,np) )
    allocate( bc_east (is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,np) )
    bc_south = 0; bc_north=0; bc_west=0; bc_east=0
    allocate( M2_psi_diss(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,np) ); M2_psi_diss = 0
  endif

  if (enable_idemix_M2) then
   allocate( E_M2(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,np,3) ); E_M2=0
   allocate(dE_M2p(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,np,3) );dE_M2p=0
   allocate( cg_M2(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx) ); cg_M2 = 0.0; 
   allocate( kdot_x_M2(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx) ); kdot_x_M2 = 0.0; 
   allocate( kdot_y_M2(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx) ); kdot_y_M2 = 0.0; 
   allocate(  forc_M2(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,np) ); forc_M2 = 0
   allocate( u_M2(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,np) ); u_M2=0.
   allocate( v_M2(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,np) ); v_M2=0.
   allocate( w_M2(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,np) ); w_M2=0.
   allocate( E_struct_M2(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); E_struct_M2 = 0
   allocate( E_M2_int(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx) ); E_M2_int = 0
  endif


  if (enable_idemix_niw) then
    allocate( omega_niw(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx) ); omega_niw=0
    allocate( E_niw(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,np,3) ); E_niw=0
    allocate(dE_niwp(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,np,3) );dE_niwp=0
    allocate( cg_niw(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx) ); cg_niw = 0.0; 
    allocate( kdot_x_niw(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx) ); kdot_x_niw = 0.0; 
    allocate( kdot_y_niw(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx) ); kdot_y_niw = 0.0; 
    allocate(  forc_niw(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,np) ); forc_niw = 0
    allocate( u_niw(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,np) ); u_niw=0.
    allocate( v_niw(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,np) ); v_niw=0.
    allocate( w_niw(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,np) ); w_niw=0.
    allocate( E_struct_niw(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); E_struct_niw = 0
    allocate( E_niw_int(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx) ); E_niw_int = 0
  endif




end subroutine allocate_idemix_module
