!=======================================================================
!
! slope convection
!
!=======================================================================

module config_module
 real*8 :: fac=2.0,visc=4e-3,mix=2e-4,L0,H0
 real*8 :: pipe_diffusivity = 5e-3
 integer :: end_of_pipe, start_of_pipe, sucked_nz, pipe_position
 real*8 :: Tp, Tm, Sp, Sm
end module config_module

subroutine set_parameter
 use main_module   ! check in directory for_src/main 
 use config_module   
 use diagnostics_module    
 implicit none 
 real*8, external :: get_drhodT, get_drhodS
 
  nx=int(32*fac)
  nz=int(32*fac)
  ny=1
! Set here time step
! cfl dt < dx/u
  dt_tracer=0.25/fac 
  dt_mom   =0.25/fac 
     
  enable_conserve_energy = .false.
  coord_degree           = .false.
! need hydrostatic approximation? Then set to true
  enable_hydrostatic     = .false.
     
  eq_of_state_type       =  1

! parameters for surface pressure integration
  congr_epsilon = 1e-6
  congr_max_iterations = 5000
  congr_epsilon_non_hydro=   1e-6
  congr_max_itts_non_hydro = 5000    

! Set here viscosities
  enable_explicit_vert_friction = .true.
  kappaM_0 = visc/fac**2 ! vertical viscosity
  enable_hor_friction = .true.
  A_H = visc/fac**2 ! lateral viscosity

! set here diffusivities
  !enable_hor_diffusion = .true.
  !kappah_0 = mix/fac**2
  !K_h = mix/fac**2

! tracer advection scheme
  enable_superbee_advection = .true.
  
  enable_tempsalt_sources = .true.

! Set here run length in seconds and snapshot (output) interval snapint     
  runlen =  1000000.0
  enable_diag_ts_monitor = .true.; ts_monint =dt_tracer
  enable_diag_snapshots  = .true.; snapint  = 5
  
  
 ! Linear equation of state: Delta rho = rho0*(betaT Delta T - betaS Delta S)
 Tp = 10. + .5/(-get_drhodT(34d0,10d0,0d0) )
 Tm = 10.
 Sp = 34. + .5/get_drhodS(34d0,10d0,0d0)
 Sm = 34.
 pipe_position = nx/3
 sucked_nz = 10 ! Grid points that T-, S- get sucked into area of T+, S+
 start_of_pipe = 20 ! Update this according to adv_main.f90
 end_of_pipe = nz - 20
   
end subroutine set_parameter


subroutine set_grid
 use main_module   
 use config_module   
 implicit none
! Set here the resolution of your model
 dxt(:)=0.5/fac ! in m
 dyt(:)=0.5/fac 
 dzt(:)=0.5/fac 
 L0 = dxt(is_pe)*nx
 H0 = dzt(1)*nz
end subroutine set_grid

subroutine set_coriolis
end subroutine set_coriolis

subroutine set_initial_conditions
 use main_module   
 use config_module   
 implicit none
 integer :: i,j, k
 
  if (pipe_position>=is_pe .and. pipe_position<=ie_pe) then 
   ! modify land mask
   maskT(pipe_position,:,start_of_pipe:end_of_pipe)=0.
  endif
 
 ! recalculate related masks 
  call border_exchg_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,maskT) 
  call setcyclic_xyz   (is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,maskT)
  maskU=maskT
  do i=is_pe-onx,ie_pe+onx-1
     maskU(i,:,:)=min(maskT(i,:,:),maskT(i+1,:,:))
  enddo
  call border_exchg_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,maskU) 
  call setcyclic_xyz   (is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,maskU)
  maskV=maskT
  do j=js_pe-onx,je_pe+onx-1
      maskV(:,j,:)=min(maskT(:,j,:),maskT(:,j+1,:))
  enddo
  call border_exchg_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,maskV) 
  call setcyclic_xyz   (is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,maskV)
  maskZ=maskT
  do j=js_pe-onx,je_pe+onx-1
   do i=is_pe-onx,ie_pe+onx-1
     maskZ(i,j,:)=min(maskT(i,j,:),maskT(i,j+1,:),maskT(i+1,j,:))
   enddo
  enddo
  call border_exchg_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,maskZ) 
  call setcyclic_xyz   (is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,maskZ)
  maskW=maskT
  do k=1,nz-1
    maskW(:,:,k)=min(maskT(:,:,k),maskT(:,:,k+1))
  enddo


 do i = 1, pipe_position ! set temp, salt left of pipe boundary
 
  if (i>=is_pe .and. i<=ie_pe) then 
    temp(i,:,:,tau) = Tp ! temp, salt above pipe
    salt(i,:,:,tau) = Sp
    do k = 1, start_of_pipe + sucked_nz
       temp(i,:,k,tau) = Tm
       salt(i,:,k,tau) = Sm
    end do
   do k = (start_of_pipe + sucked_nz + 1), end_of_pipe
      temp(i,:,k,tau) = Tp
      salt(i,:,k,tau) = Sp
   end do
  endif
  
 end do
   
 do i= (pipe_position+1), nx ! set temp, salt to the right  of the pipe
  if (i>=is_pe .and. i<=ie_pe) then
    temp(i,:,:,tau) = Tp
    salt(i,:,:,tau) = Sp
     do k = 1, start_of_pipe
        temp(i,:,k,tau) = Tm ! Set temp, salt at bottom half
        salt(i,:,k,tau) = Sm
     end do
   endif  
 end do

 call border_exchg_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,temp(:,:,:,tau)) 
 call setcyclic_xyz   (is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,temp(:,:,:,tau))
 call border_exchg_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,salt(:,:,:,tau)) 
 call setcyclic_xyz   (is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,salt(:,:,:,tau))
 
 temp(:,:,:,taum1) = temp(:,:,:,tau)
 salt(:,:,:,taum1) = salt(:,:,:,tau)
end subroutine set_initial_conditions


subroutine set_forcing
 use main_module   
 use config_module   
 implicit none
 real*8 :: flx(js_pe-onx:je_pe+onx,nz)
 
 ! diffusive flux across pipe boundaries
 ! does not work of pipe_position is a PE boundaries
 if (pipe_position==is_pe .or. pipe_position==ie_pe) then
  if (my_pe==0) print*,' ERROR: pipe boundary is at PE boundary '
  if (my_pe==0) print*,'        change pipe position or number of PEs '
  call halt_stop(' in set_forcing')
 endif
 
 if (pipe_position>=is_pe .and. pipe_position<=ie_pe) then ! also for parallel use
   flx = pipe_diffusivity*(temp(pipe_position+1,:,:,taum1) - temp(pipe_position-1,:,:,taum1) )/(2*dxt(pipe_position)) 
   temp_source(pipe_position-1,:,:) = +flx
   temp_source(pipe_position+1,:,:) = -flx  
 endif
 call border_exchg_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,temp_source) 
 call setcyclic_xyz   (is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,temp_source)
end subroutine set_forcing


subroutine set_topography
 use main_module   
 use config_module   
 implicit none

! kbot is array controlling topography
! kbot=1 everywhere means no topographic features
! kbot=0 means land point
 kbot =1
end subroutine set_topography


subroutine set_diagnostics
end subroutine set_diagnostics


subroutine set_particles
end subroutine set_particles

