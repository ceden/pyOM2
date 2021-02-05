
!=======================================================================
! 1 D test case for tke model
!======================================================================= 

module config_module
 ! use this module only locally in this file
 implicit none
 real*8,parameter :: N_0     = 0.02
end module config_module


subroutine set_parameter
 ! ----------------------------------
 !       set here main parameter
 ! ----------------------------------
 use main_module   
 use config_module
 use diagnostics_module   
 use tke_module   

 implicit none
 nx   = 1; nz   = 50; ny   = 1
 dt_mom     = 360
 dt_tracer  = 360

 coord_degree           = .false.
 enable_cyclic_x        = .true.
 enable_cyclic_y        = .true.
 enable_hydrostatic     = .true.
 eq_of_state_type       = 5

 congr_epsilon = 1e-12
 congr_max_iterations = 5000
 !enable_streamfunction = .true.
     
 runlen =  365*86400.*10
 enable_diag_ts_monitor = .true.; ts_monint = dt_tracer
 enable_diag_snapshots  = .true.; snapint  =  86400./6.
 !enable_diag_energy     = .true.; energint =  snapint; energfreq = dt_tracer

 enable_implicit_vert_friction = .true.;
 enable_tke = .true.
 alpha_tke = 30.0
 tke_mxl_choice = 2
 
end subroutine set_parameter



subroutine set_grid
 use main_module   
 use config_module   
 implicit none
 dxt    = 5e3
 dyt    = 5e3
 dzt    = 1.0
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
 use config_module
 implicit none
 integer :: i,j,k
 real*8 :: alpha,get_drhodT

 salt = 35.0
 alpha = get_drhodT(35d0,0d0,0d0)
 do k=1,nz
   do j=js_pe-onx,je_pe+onx
     do i=is_pe-onx,ie_pe+onx
       temp(i,j,k,:)  = ( 32+rho_0/grav/alpha*(-N_0**2*zt(k))  )*maskT(i,j,k)
       !u(i,j,k,:)     = 0.5*tanh( (zt(k)-zt(nz/2))/zt(1)*3)
     enddo
   enddo
 enddo

 surface_taux = .1e-3*maskU(:,:,nz)
 surface_tauy = 0.0
 
 if (enable_tke ) then
  tke=1d-12
  do j=js_pe-onx+1,je_pe+onx
   do i=is_pe-onx+1,ie_pe+onx
     forc_tke_surface(i,j) = sqrt( (0.5*(surface_taux(i,j)+surface_taux(i-1,j)))**2  &
                                  +(0.5*(surface_tauy(i,j)+surface_tauy(i,j-1)))**2 )**(3./2.)
   enddo
  enddo
 endif

end subroutine set_initial_conditions


subroutine set_forcing
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





