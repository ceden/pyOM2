
!=======================================================================
!  idealised Southern Ocean, same as in Viebahn and Eden (2010) Ocean modeling
!======================================================================= 


module config_module
 ! use this module only locally in this file
 implicit none
 real*8,parameter :: hRESOLVE = 0.5 ! 1 in original model
 real*8 :: L_y,L_x
 real*8 :: t_rest=30*86400
 real*8:: phi0 = 30.0 /180. *3.1415
end module config_module


subroutine set_parameter
 ! ----------------------------------
 !       set here main parameter
 ! ----------------------------------
 use main_module   
 use config_module
 use diagnostics_module   

 implicit none
 nx   = int( 64*hRESOLVE ); nz   = 3; ny   = int( 64*hRESOLVE )
 dt_mom     = 200.0/hRESOLVE
 dt_tracer  = 200.0/hRESOLVE !*5

 coord_degree           = .false.
 enable_cyclic_x        = .true.
 enable_hydrostatic     = .true.
 eq_of_state_type       = 1

 congr_epsilon = 1e-6
 congr_max_iterations = 5000
 enable_streamfunction = .true.
     
 !enable_upwind3_advection = .true.
 enable_dst3_advection = .true.
 !enable_superbee_advection = .true.
 enable_AB_time_stepping   = .false.
 enable_momentum_equation  = .false.

 !kappah_0=1.e-4/vRESOLVE
 !kappam_0=1.e-3/vRESOLVE
 enable_conserve_energy = .false.
 runlen =  365*86400.*0.5
 enable_diag_ts_monitor = .true.; ts_monint = dt_tracer
 enable_diag_snapshots  = .true.; snapint  =  86400./20.
end subroutine set_parameter



subroutine set_grid
 use main_module   
 use config_module   
 implicit none
 dxt    = 5e3/hRESOLVE
 dyt    = 5e3/hRESOLVE
 dzt    = 100.0!/vRESOLVE
end subroutine set_grid

subroutine set_coriolis
end subroutine set_coriolis


subroutine set_initial_conditions
 use main_module
 use tke_module
 use idemix_module
 use config_module
 use obc_module
 implicit none
 integer :: i,j,k
 real*8 :: alpha,get_drhodT,L_R,k0

 salt = 35.0
 temp = 0.0
 u(:,:,:,1)=10*maskU
 u(:,:,:,2)=10*maskU
 u(:,:,:,3)=10*maskU

 alpha = get_drhodT(35d0,0d0,0d0) 
 L_R = 25e3!(nz*dzt(1)*N_0)/(pi*coriolis_t(i,j) )
 k0 = 2*2*pi/(nx*dxt(is_pe))
 do k=1,nz
  do j=js_pe-onx,je_pe+onx
   do i=is_pe-onx,ie_pe+onx
    !temp(i,j,k,:)    =exp(-(xt(i)-L_x/2.)**2/L_R**2-(yt(j)-L_y/2.)**2/L_R**2)*maskT(i,j,k)
    temp(i,j,k,:)    =sin(xt(i)*k0)*maskT(i,j,k)
   enddo
  enddo
 enddo

end subroutine set_initial_conditions



subroutine set_forcing
end subroutine set_forcing

subroutine set_topography
 use main_module   
 use config_module   
 implicit none

 L_y = 0.0; if (my_blk_j == n_pes_j) L_y = yu(ny)
 call global_max(L_y)
 L_x = 0.0; if (my_blk_i == n_pes_i) L_x = xu(nx)
 call global_max(L_x)
 if (my_pe==0) print*,' domain size is ',L_x,' m x ',L_y,' m'

 kbot=1
end subroutine set_topography




subroutine set_diagnostics
end subroutine set_diagnostics

subroutine set_particles
end subroutine set_particles





