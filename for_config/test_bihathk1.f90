
!=======================================================================
!  test for isopycnal mixing
!======================================================================= 
module config_module
 ! use this module only locally in this file
 implicit none
   real*8 :: Ly=1000e3,Lz = 2000. ,f0 = 1e-4 ,N0=30e-4
end module config_module 


subroutine set_parameter
 ! ----------------------------------
 !       set here main parameter
 ! ----------------------------------
 use main_module     
 use isoneutral_module  
 use biharmonic_thickness_module
 use diagnostics_module 
 use config_module 
 implicit none

  nx   = 1; nz   = 40; ny  = 30
  dt_mom    = 1200
  dt_tracer = dt_mom

 coord_degree           =  .false.
 eq_of_state_type       =  5
 enable_hydrostatic      = .true.
 enable_cyclic_x         = .true.

 congr_epsilon = 1e-8
 congr_max_iterations = 5000
 
  runlen = dt_mom*1000000 !365*86400.*2000  

  enable_diag_ts_monitor = .true.; ts_monint = dt_mom!/12.
  enable_diag_snapshots  = .true.; snapint  =  86400.!*30
  !enable_diag_tracer_content = .true.; trac_cont_int=dt_mom
  
  
 enable_conserve_energy             = .true.
 enable_store_lateral_friction_heat = .true.

 enable_superbee_advection = .true.


  !enable_neutral_diffusion = .true.; 
      K_iso_0 = 00.0
      K_iso_steep = 00.0
      iso_dslope=0.002
      iso_slopec=0.01
  K_gm_0 = 2000.0    
  enable_skew_diffusion = .true.

 enable_biharmonic_friction    =  .true.
 A_hbi    =  1e12
 
 
 enable_biharmonic_thickness_bolus_skew = .true.
 A_thkbi = A_hbi
 mld0_thk = 0.
 biha_dslope=1d-12
 biha_slopec=4d-12 
  
 !enable_biharmonic_thickness_backscatter_skew = .true.
 !biharmonic_thickness_backscatter_smooth = 1
 !A_thk_0 = -200.0
 !enable_biharmonic_thickness_backscatter_integrate_energy = .true.

 !enable_biharmonic_thickness_explicit = .true.
 !A_thkbi = A_hbi
 !biharmonic_thickness_mixing_iter = 30
 
end subroutine set_parameter



subroutine set_grid
  use main_module   
  use config_module   
  implicit none
   dyt = Ly/ny
   dxt = dyt(js_pe)
   dzt = Lz/nz 
   
end subroutine set_grid


subroutine set_coriolis
 use main_module   
 use config_module   
 implicit none
  coriolis_t = f0
end subroutine set_coriolis


subroutine set_initial_conditions
 use main_module   
 use config_module    
 implicit none
 integer :: i,j,k
 real*8 :: y,alpha,get_drhodT
 

 alpha = get_drhodT(35d0,5d0,0d0)  
 do k=1,nz
  do j=js_pe-onx,je_pe+onx
   do i=is_pe-onx,ie_pe+onx
    y = yt(j)/Ly
    salt(i,j,k,:) = 35
    temp(i,j,k,:) = 20-N0**2*zt(k)/grav/alpha*rho_0 +  2*(tanh((y-0.5)*20)+1.0)
   enddo
  enddo
 enddo
 temp(:,:,nz,:) = temp(:,:,nz-3,:)
 temp(:,:,nz-1,:) = temp(:,:,nz-3,:)
 temp(:,:,nz-2,:) = temp(:,:,nz-3,:)

end subroutine set_initial_conditions

subroutine set_momentum_forcing
end subroutine set_momentum_forcing

subroutine set_forcing
end subroutine set_forcing

subroutine set_topography
 use main_module   
 implicit none
 !integer :: j
 kbot=1
 
 !do j=js_pe-onx,je_pe+onx
 !  if (j>=ny/4.and.j<=3*ny/4)   kbot(:,j) = nz/5
 !enddo
end subroutine set_topography

subroutine set_diagnostics
end subroutine set_diagnostics

subroutine set_particles
end subroutine set_particles