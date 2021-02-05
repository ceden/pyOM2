
!=======================================================================
! Kelvin Helmholtz instability
!=======================================================================



subroutine set_parameter
 use main_module   
 use diagnostics_module
 implicit none
 
  nx=64;ny=1;nz=nx
 
  dt_mom=0.1
  dt_tracer = dt_mom

  enable_conserve_energy = .false. 
  coord_degree         =.false.
  enable_cyclic_x      =.true.
  enable_cyclic_y      =.true.
  enable_hydrostatic   =.false.
  eq_of_state_type = 1 

  congr_epsilon = 1d-6
  congr_epsilon_non_hydro=   1d-4
  

  enable_explicit_vert_friction = .true.;
  enable_hor_friction = .true.;
  enable_hor_diffusion = .true.; 
             
  A_h = 1e-3
  K_h = A_h
  kappaM_0 = A_h
  kappaH_0 = A_h
         
  

  
  enable_diag_snapshots  = .true.; snapint  = 5.0
  enable_diag_ts_monitor = .true.; ts_monint = snapint
  snapint = 5.0 !dt*10
  runlen = snapint*400

end subroutine set_parameter


subroutine set_grid
 use main_module   

 implicit none
 dxt(:)=0.25
 dyt(:)=0.25
 dzt(:)=0.25
end subroutine set_grid

subroutine set_coriolis
end subroutine set_coriolis





subroutine set_initial_conditions
  use main_module   
  implicit none
  integer :: i,j,k
  real*8 :: Lx,z0,z1
  
  Lx = nx*dxt(is_pe)
  z0 =  zt(nz/2)
  z1 = 0.25/2.
  do k=1,nz
   do j=js_pe,je_pe
    do i=is_pe-onx,ie_pe+onx  
     u(i,j,k,tau)    = 0.6+0.5*tanh( (zt(k)-z0)/z1) 
     temp(i,j,k,tau) = 10*(1+tanh( (zt(k)-z0)/z1) +0.2*sin( xt(i)/Lx*12*pi) )              
    enddo
   enddo
  enddo  
    
end subroutine set_initial_conditions




subroutine set_forcing

end subroutine set_forcing



subroutine set_topography
 use main_module   
 implicit none
 kbot =1
end subroutine set_topography


subroutine set_diagnostics
end subroutine set_diagnostics


subroutine set_particles
end subroutine set_particles
