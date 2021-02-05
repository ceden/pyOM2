
!=======================================================================
! 1 D test case for idemix and tke model
!======================================================================= 

module config_module
 ! use this module only locally in this file
 implicit none
 real*8,parameter :: N_0 = 0.008, z_b = 800.0  ! N=N_0 exp(z/z_b)
 real*8, parameter :: U00 = 0.5,V00 = 0.0, zu0=200.0,zu1=500, zu2 = 1200
 real*8, allocatable :: u_ini(:,:,:),v_ini(:,:,:)
end module config_module


subroutine set_parameter
 ! ----------------------------------
 !       set here main parameter
 ! ----------------------------------
 use main_module   
 use config_module
 use diagnostics_module   
 use tke_module   
 use idemix_module   

 implicit none
 !nx   = 1; nz   = 200; ny   = 1
 nx   = 1; nz   = 50; ny   = 1
 dt_mom     = 360
 dt_tracer  = dt_mom

 coord_degree           = .false.
 enable_cyclic_x        = .true.
 enable_cyclic_y        = .true.
 enable_hydrostatic     = .true.
 eq_of_state_type       = 5 ! 5????

 !enable_momentum_equation = .false.
 enable_thermodynamic_equation = .false.
 
 congr_epsilon = 1e-12
 congr_max_iterations = 5000
 !enable_streamfunction = .true.
     
 runlen =  365*86400.*10
 enable_diag_ts_monitor = .true.; ts_monint = dt_tracer
 enable_diag_snapshots  = .true.; snapint  =  3600.
 enable_diag_energy     = .true.; energint =  snapint; energfreq = dt_tracer
 
 !enable_conserve_energy = .false.
 enable_store_cabbeling_heat = .true.
 enable_store_bottom_friction_tke = .true.

 enable_momentum_sources_with_AB = .true.

 enable_implicit_vert_friction = .true.;
 enable_tke = .true.
 
 enable_idemix = .true.
 enable_idemix3 = .true.
 tau_v = 3*86400.
 mu0 = 4./3.
 mu0_min = 0.
 tc_max = 1e-3
 Noverf_min0 = 0.
 !AB_eps = 0.03
 
 idemix3_leewaves_intervall = dt_tracer
  enable_leewaves = .true.
  nk  = 32
  nph = 32
  
end subroutine set_parameter



subroutine set_grid
 use main_module   
 use config_module   
 implicit none
 dxt    = 5e3
 dyt    = 5e3
 dzt    = 2000./nz !10.0
end subroutine set_grid

subroutine set_coriolis
 use main_module   
 use config_module
 implicit none
 coriolis_t = 2*omega*sin(30.0 /180.*pi)
end subroutine set_coriolis


subroutine set_initial_conditions
 use main_module
 use tke_module
 use idemix_module   
 use config_module
 implicit none
 integer :: i,j,k
 real*8 :: alpha,get_drhodT,N2


 allocate( u_ini(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); u_ini=0.
 allocate( v_ini(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); v_ini=0.

 salt = 35.0
 temp(:,:,nz,:)=25d0
 do k=nz-1,1,-1
   do j=js_pe-onx,je_pe+onx
      do i=is_pe-onx,ie_pe+onx  ! b_z=N^2 , b=-rho/rho_0 g,  T=rho/rho_T, T=-b rho_0/( rho_T g)   
       alpha = get_drhodT(salt(i,j,k,1),temp(i,j,k+1,1),abs(zt(k)))
       N2 = ( N_0*exp(zt(k)/z_b)  )**2
       temp(i,j,k,:)  = ( temp(i,j,k+1,:) +dzw(k)*rho_0/(grav*alpha)*N2 )*maskT(i,j,k)
     enddo
   enddo
 enddo

 do k=1,nz
   u_ini(:,:,k)     = 0.1*(1-tanh( (zt(k)+zu1)/zu0))
   !v_ini(:,:,k)     = V00*tanh( (zt(k)+zu1)/zu0)
   !u_ini(:,:,k)     = 0.1+U00*exp( -(zt(k)+zu1)**2/zu0**2)  
   !v_ini(:,:,k)     = V00*exp( -(zt(k)+zu1)**2/zu0**2)  
 enddo
 
 u(:,:,:,tau)   = u_ini
 v(:,:,:,tau)   = v_ini
 u(:,:,:,taum1) = u_ini
 v(:,:,:,taum1) = v_ini
 u(:,:,:,taup1) = u_ini
 v(:,:,:,taup1) = v_ini

 if (enable_idemix ) then
     !do k=1,nz
     !  E_iw(:,:,k,:)     = 1e-3*exp( -(zt(k)-zt(nz/2))**2/200.0**2)
     !enddo
  if (enable_idemix3 ) then
     ! E+ = 0.5(E_s+E_d)
     ! E- = 0.5(E_s-E_d)

     do k=1,nz
    !   F_s(:,:,k,:)     = 1e-3*exp( -(zt(k)+1200)**2/200.0**2)
    !   G_s(:,:,k,:)     = .1e-3*exp( -(zt(k)+1200)**2/200.0**2)
     enddo
     !forc_iw_bottom_u =  0.5e-6
     !forc_iw_bottom_v =  0.5e-6
     !forc_iw_surface_u =  0.5e-6
     !forc_iw_surface_v =  0.5e-6
     !forc_iw_surface = 2*.5e-6
     
    if (enable_leewaves) then
      h_rms = 50.
      ph_s = 0
      k_s = 2*sqrt(2*(nu+0.5))/20e3
      k_n = k_s
    else
      forc_iw_bottom_u =  1.5e-5
      forc_iw_bottom_v =  1.5e-5  
    endif       
     
  else
     forc_iw_bottom =  2*0.5e-6
  endif
endif

end subroutine set_initial_conditions


subroutine set_forcing
 ! u_t =   f v - f v*  + Fu_z
 ! v_t = - f u + f u*  + Fv_z
 use main_module   
 use config_module
 implicit none
 integer :: k
 do k=1,nz
    u_source(:,:,k)     = -coriolis_t(:,:)*v_ini(:,:,k)
    v_source(:,:,k)     =  coriolis_t(:,:)*u_ini(:,:,k)
 enddo
end subroutine set_forcing

subroutine set_topography
 use main_module   
 implicit none
 kbot=1
end subroutine set_topography


subroutine set_diagnostics
end subroutine set_diagnostics

subroutine set_particles
end subroutine set_particles





