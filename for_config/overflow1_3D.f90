

!=======================================================================
! 3D slope convection 
!=======================================================================

module config_module
 real*8 :: fac=2.0,mix=2e-06,L0,H0,B0
 real*8, allocatable :: T0(:,:,:)
end module config_module

subroutine set_parameter
 use main_module   
 use config_module   
 use diagnostics_module  
 use tke_module   
 implicit none
  nx=int(60*fac)
  nz=int(20*fac)
  ny=int(60*fac)

  dt_tracer=0.25/fac
  dt_mom   =0.25/fac
     
  enable_conserve_energy = .false.
  coord_degree           = .false.
  enable_hydrostatic     = .false.
     
  eq_of_state_type       =  1

  congr_epsilon = 1e-8
  congr_max_iterations = 5000
  congr_epsilon_non_hydro=   1e-8
  congr_max_itts_non_hydro = 5000

  enable_bottom_friction = .true.
  r_bot = 0.7

  enable_explicit_vert_friction = .true.
  kappam_0 = mix/fac**2
  !vertikale reibung
  enable_hor_friction = .true.
  a_h = mix/fac**2
  !horizontale reibung

  enable_superbee_advection = .true.
  !enable_hor_diffusion = .true
  !kappah_0 = mix/fac**2
  !k_h = mix/fac**2
     
  !enable_tempsalt_sources =  .true. ! restoring zone


  runlen =  86400
  enable_diag_ts_monitor = .true.; ts_monint =dt_tracer
  enable_diag_snapshots  = .true.; snapint  = 1
end subroutine set_parameter


subroutine set_grid
 use main_module   
 use config_module   
 implicit none
 dxt(:)=0.01/fac
 dyt(:)=0.01/fac
 dzt(:)=0.01/fac
 
 L0 = dxt(is_pe)*nx !länge tank
 H0 = dzt(1)*nz     !höhe tank
 B0 = dyt(js_pe)*ny !breite tank
end subroutine set_grid

subroutine set_coriolis
use main_module   
implicit none
!coriolis_t = 0.1047 !1rpm
!coriolis_t = 0.2094 !2rpm
coriolis_t = 0.4189 !4rpm
!coriolis_t = 0.8378 !8rpm
end subroutine set_coriolis



subroutine set_initial_conditions
 use main_module   
 use config_module   
 implicit none
 integer :: i,slope
 real*8 :: temp_shelf,shelf

 allocate( T0(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); T0=0.0
 slope = nx/4 !Wo soll der Hang beginnen? Muss Integer sein!
 shelf = H0/4 !Wie tief soll das Schelf sein?
 temp_shelf = 10
 T0 = 12 
 do i=is_pe-onx,ie_pe+onx
    if (i<=slope) T0(i,:,:)=temp_shelf
 enddo
 temp(:,:,:,taum1) = T0
 temp(:,:,:,tau)   = T0
end subroutine set_initial_conditions


subroutine set_forcing
 use main_module   
 use config_module   
 implicit none
 integer :: i
 !if (enable_tempsalt_sources) then
 ! do i=1,nx/4
 !  if (i>=is_pe.and.i<=ie_pe) temp_source(i,:,:)=maskT(i,:,:)*(T0(i,:,:)-temp(i,:,:,taum1))/(dt_mom*10)
 ! enddo
 !endif
end subroutine set_forcing




subroutine set_topography
 use main_module   
 use config_module   
 implicit none
 integer :: i,k,j,slope
 real*8 :: fxa,shelf,gap1,gap2,a,b
 kbot =1

 slope = nx/4 !Wo soll der Hang beginnen? Muss Integer sein!
 shelf = H0/4 !Wie tief soll das Schelf sein?
 gap1 = nx*2/3 !Wo soll die Lücke in der Barriere beginnen?
 gap2 = nx*5/6 !Wo soll sie enden?
 b=(nx*shelf-H0*slope)/(nx-slope)
 a=(shelf-b)/slope
 
 do i=is_pe,ie_pe 
     if (i<=slope) fxa=shelf
     if (i>slope)  fxa=a*i+b
     do k=1,nz
       if (zt(k).lt.-fxa) kbot(i,:)=k
     enddo
     if (i == slope) then
          do j=js_pe,je_pe 
             if (j<gap1) then
                kbot(i,j)=0
             endif
             if (j>gap2) then
                kbot(i,j)=0
             endif
          enddo
     end if
 enddo
end subroutine set_topography

subroutine set_diagnostics
end subroutine set_diagnostics

subroutine set_particles
end subroutine set_particles



