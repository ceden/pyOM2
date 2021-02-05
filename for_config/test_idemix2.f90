
!=======================================================================
! 2 D test case for idemix and tke model
!======================================================================= 

module config_module
 ! use this module only locally in this file
 implicit none
 real*8,parameter :: N_0 = 0.008, z_b = 800.0  ! N=N_0 exp(z/z_b)
 real*8, parameter :: U00 = 0.5, y_star=50.e3,dy_star=20e3, h_star=500, dh_star = 200
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
 nx   = 1; nz   = 50; ny   = 20
 dt_mom     = 1800
 dt_tracer  = dt_mom

 coord_degree           = .false.
 enable_cyclic_x        = .true.
 !enable_cyclic_y        = .true.
 enable_hydrostatic     = .true.
 eq_of_state_type       = 1 

 !enable_momentum_equation = .false.
 enable_thermodynamic_equation = .false.
 
 congr_epsilon = 1e-12
 congr_max_iterations = 5000
 !enable_streamfunction = .true.
     
 runlen =  365*86400.*10
 enable_diag_ts_monitor = .true.; ts_monint = dt_mom
 enable_diag_snapshots  = .true.; snapint  =  dt_mom
 !enable_diag_energy     = .true.; energint =  snapint; energfreq = dt_tracer
 
 !enable_conserve_energy = .false.
 enable_store_cabbeling_heat = .true.
 enable_store_bottom_friction_tke = .true.

 enable_momentum_sources_with_AB = .true.

 enable_hor_friction = .true.
 A_h = 200

 !enable_implicit_vert_friction = .true.;
 !enable_tke = .true.
 
 enable_idemix = .true.
 enable_idemix3 = .true.
 tau_v = 3*86400.
 mu0 = 4./3.
 mu0_min = 0.
 tc_max = 1e-3
 Noverf_min0 = 0.
 !AB_eps = 0.03
end subroutine set_parameter



subroutine set_grid
 use main_module   
 use config_module   
 implicit none
 dxt    = 5e3
 dyt    = 100e3/ny 
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
 real*8 :: aloc(nx,ny),bloc(ny,nz)
 real*8 :: uz(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz)

 alpha = get_drhodT(35d0,5d0,0d0) 

 allocate( u_ini(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); u_ini=0.
 allocate( v_ini(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); v_ini=0.
 salt = 35.0
 temp(:,:,nz,:)=25d0
 do k=nz-1,1,-1
   do j=js_pe-onx,je_pe+onx
      do i=is_pe-onx,ie_pe+onx  ! b_z=N^2 , b=-rho/rho_0 g,  T=rho/rho_T, T=-b rho_0/( rho_T g)   
       N2 = ( N_0*exp(zt(k)/z_b)  )**2
       temp(i,j,k,:)  = ( temp(i,j,k+1,:) +dzw(k)*rho_0/(grav*alpha)*N2 )*maskT(i,j,k)
     enddo
   enddo
 enddo

 do j=js_pe,je_pe
  do k=1,nz
   u_ini(:,j,k)     = U00*exp(-(yt(j)-y_star)**2/dy_star**2)*0.5*(1+tanh( (zt(k)+h_star)/dh_star ))
   !uz   (:,j,k)     = U00/dh_star*exp(-(yt(j)-y_star)**2/dy_star**2)*(1-tanh( (zt(k)+h_star)/dh_star )**2)
  enddo
 enddo
 do k=2,nz-1
   uz(:,:,k) = (u_ini(:,:,k+1)-u_ini(:,:,k-1))/(dzt(k+1)+dzt(k))
 enddo
 uz(:,:,1)=uz(:,:,2)
 uz(:,:,nz)=uz(:,:,nz-1)


 u(:,:,:,tau)   = u_ini; !v(:,:,:,tau)   = v_ini
 u(:,:,:,taum1) = u_ini; !v(:,:,:,taum1) = v_ini
 u(:,:,:,taup1) = u_ini; !v(:,:,:,taup1) = v_ini

 ! rho0 fu = -p_y, p_z = -g rho, rho0 f u_z = -g rho_y,  rho_y = - rho0 f u_z/g = alpha T_y
 aloc(1,js_pe:je_pe) = dyt(js_pe:je_pe) 
 call pe0_recv_2D(nx,ny,aloc)
 bloc(1,:)=0.
 do j=1,ny-1
    bloc(j+1,:)=bloc(j,:) + aloc(1,j)*uz(1,j,:)*coriolis_t(is_pe,j)/grav/alpha*rho_0
    call pe0_bcast(bloc(j+1,:),nz)
 enddo
 do j=js_pe,je_pe
   do k=1,nz
    temp(:,j,k,:)  = temp(:,j,k,:) + bloc(j,k)
   enddo
 enddo


 if (enable_idemix.and.enable_idemix3 ) then
    forc_iw_surface_u =  0.5e-6
    forc_iw_surface_v =  0.5e-6
 endif

end subroutine set_initial_conditions


subroutine set_forcing
 ! u_t =   f v - f v* -p_x + Fu_z
 ! v_t = - f u + f u* -p_y + Fv_z
 use main_module   
 use config_module
 implicit none
 integer :: k
 do k=1,nz
    !u_source(:,:,k)     = -coriolis_t(:,:)*v_ini(:,:,k)
    !v_source(:,:,k)     =  coriolis_t(:,:)*u_ini(:,:,k)
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





