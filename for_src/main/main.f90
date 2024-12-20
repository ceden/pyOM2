




program main
!=======================================================================
!      Top level driver 
!=======================================================================
      use main_module   
      use eke_module   
      use tke_module   
      use idemix_module   
      use diagnostics_module   
      use timing_module   
      implicit none
      integer :: ierr,otaum1,iargc,n
      character (len=80) :: arg
!---------------------------------------------------------------------------------
!       Initialize MPI setting
!---------------------------------------------------------------------------------
      call mpi_init(ierr)
      call my_mpi_init(my_comm)
      if (my_pe==0) print '(/,a,f5.3)' ,'here is pyOM2 version ',version

      n_pes_i = 1; n_pes_j = 1
      if (n_pes>1) then
       if (iargc() < 2) then
        call halt_stop(' not enough command line input')
       endif
       call getarg(1,arg); read(arg,*) n_pes_i
       call getarg(2,arg); read(arg,*) n_pes_j
       if (my_pe==0) print'(/a,i4,a,i4,a)','using ',n_pes_i,' x ',n_pes_j ,' PEs'
      endif
      call tic('setup')
!---------------------------------------------------------------------------------
!      Initialize model
!---------------------------------------------------------------------------------
      itt=0
      call setup
!---------------------------------------------------------------------------------
!     read restart if present
!---------------------------------------------------------------------------------
      if (my_pe==0) print'(/,a)','Reading restarts:'
      call read_restart(itt)
      call read_tracer_restart
      if (enable_diag_averages)     call diag_averages_read_restart
      if (enable_diag_variances)    call diag_variance_read_restart
      if (enable_diag_energy)       call diag_energy_read_restart
      if (enable_diag_overturning)  call diag_over_read_restart
      if (enable_diag_particles)    call particles_read_restart
!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
      enditt = itt+int(runlen/dt_tracer)
      if (my_pe == 0 ) then
        print'(/,a,e8.2,a)','Starting integration for ',runlen,' s'
        print'(a,i10,a,i10,/)',' from time step ',itt,' to ',enditt
      endif
      call toc('setup')
!---------------------------------------------------------------------------------
!      Begin main model loop
!---------------------------------------------------------------------------------
      
      do while (itt < enditt) 
        call tic('main loop')
        call set_forcing

        !if (enable_idemix) call set_idemix_parameter
        !if (enable_idemix_M2 .or. enable_idemix_niw) call set_spectral_parameter
        call tic('parameter')
        call set_idemix_parameter
        call set_eke_diffusivities
        call set_tke_diffusivities
        call toc('parameter')

        call tic('rossmix')
        call rossmix_main
        call rossmix2_main
        call toc('rossmix')

        call tic('mom')
        if (enable_momentum_equation) call momentum
        call toc('mom')

        call tic('temp')
        if (enable_thermodynamic_equation) then
           call thermodynamics
           call integrate_tracer
        endif
        call toc('temp')

        call tic('rossmix')        
        call rossmix2_integrate
        call toc('rossmix')

        if (enable_eke .or. enable_tke .or. enable_idemix) call calculate_velocity_on_wgrid

        call tic('eke')
        if (enable_eke)    call integrate_eke
        call toc('eke')

        call tic('idemix')
        call integrate_idemix
        call toc('idemix')

        call tic('tke')
        if (enable_tke)    call integrate_tke
        call toc('tke')

        !---------------------------------------------------------------------------------
        ! Main boundary exchange
        ! for density, temp and salt this is already done in thermodynamics.f90, tracer also
        !---------------------------------------------------------------------------------
        call border_exchg_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,u(:,:,:,taup1)) 
        call setcyclic_xyz   (is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,u(:,:,:,taup1))
        call border_exchg_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,v(:,:,:,taup1)) 
        call setcyclic_xyz   (is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,v(:,:,:,taup1))
        call set_obc_boundary_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,u(:,:,:,taup1)) 
        call set_obc_boundary_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,v(:,:,:,taup1)) 
        
        if (enable_tke) then
         call border_exchg_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,tke(:,:,:,taup1)) 
         call setcyclic_xyz   (is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,tke(:,:,:,taup1))
         call set_obc_boundary_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,tke(:,:,:,taup1)) 
        endif
        if (enable_eke) then
          call border_exchg_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,eke(:,:,:,taup1)) 
          call setcyclic_xyz   (is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,eke(:,:,:,taup1))
          call set_obc_boundary_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,eke(:,:,:,taup1)) 
        endif
        if (enable_idemix) then
           call border_exchg_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,E_iw(:,:,:,taup1)) 
           call setcyclic_xyz   (is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,E_iw(:,:,:,taup1))
           call set_obc_boundary_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,E_iw(:,:,:,taup1)) 
        endif
        if (enable_idemix_M2) then
           call border_exchg_xyp(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,np,E_M2(:,:,:,taup1)) 
           call setcyclic_xyp   (is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,np,E_M2(:,:,:,taup1))
           call set_obc_boundary_xyp(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,np,E_M2(:,:,:,taup1)) 
        endif
        if (enable_idemix_niw) then
           call border_exchg_xyp(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,np,E_niw(:,:,:,taup1)) 
           call setcyclic_xyp   (is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,np,E_niw(:,:,:,taup1))
           call set_obc_boundary_xyp(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,np,E_niw(:,:,:,taup1)) 
        endif
 
        ! diagnose vertical velocity at taup1
        if (enable_hydrostatic) call vertical_velocity
        call toc('main loop')
 
        call tic('diag')
        call diagnose
        call toc('diag')
 
        ! shift time 
        otaum1=taum1; taum1= tau; tau  = taup1; taup1= otaum1; itt=itt+1
!---------------------------------------------------------------------------------
!       End main model loop
!---------------------------------------------------------------------------------
      enddo
      if (my_pe==0) print'(/a)','end of integration'
      call write_restart(itt)
      call write_tracer_restart
      if (enable_diag_averages)    call diag_averages_write_restart
      if (enable_diag_variances)   call diag_variance_write_restart
      if (enable_diag_energy)      call diag_energy_write_restart
      if (enable_diag_overturning) call diag_over_write_restart
      if (enable_diag_particles  ) call particles_write_restart

      if (my_pe==0) then
       call get_free_iounit(n,ierr)
       open(n,file='ritt',form='formatted',status='unknown')
       write(n,*) itt
       close(n)
      endif
!--------------------------------------------------------------
!     show timing results here
!--------------------------------------------------------------
      do n = 0,n_pes
       call fortran_barrier
       if (my_pe == n) then
        print'(/,a,i4)','Timing summary for PE #',my_pe 
        print'(a,f12.1,a)',' costs for measuring      = ',timing_secs('tictoc'),' s'
        print'(a,f12.1,a)',' setup time summary       = ',timing_secs('setup'),' s'
        print'(a,f12.1,a)',' main loop time summary   = ',timing_secs('main loop') ,' s'
        print'(a,f12.1,a)','     setting parameter    = ',timing_secs('parameter') ,' s'
        print'(a,f12.1,a)','     momentum             = ',timing_secs('mom') ,' s'
        print'(a,f12.1,a)','       pressure           = ',timing_secs('press') ,' s'
        print'(a,f12.1,a)','       friction           = ',timing_secs('fric') ,' s'
        print'(a,f12.1,a)','     thermodynamics       = ',timing_secs('temp') ,' s'
        print'(a,f12.1,a)','       lateral mixing     = ',timing_secs('iso') ,' s'
        print'(a,f12.1,a)','       vertical mixing    = ',timing_secs('vmix') ,' s'
        print'(a,f12.1,a)','       equation of state  = ',timing_secs('eq_of_state') ,' s'
        print'(a,f12.1,a)','     EKE                  = ',timing_secs('eke') ,' s'
        print'(a,f12.1,a)','     IDEMIX               = ',timing_secs('idemix') ,' s'
        print'(a,f12.1,a)','       leewave_flux       = ',timing_secs('leewave_flux') ,' s'
        print'(a,f12.1,a)','       integrate_leewaves = ',timing_secs('integrate_leewaves') ,' s'
        print'(a,f12.1,a)','     TKE                  = ',timing_secs('tke') ,' s'
        print'(a,f12.1,a)','     ROSSMIX              = ',timing_secs('rossmix') ,' s'
        print'(a,f12.1,a)',' diagnostics              = ',timing_secs('diag') ,' s'
       endif
      enddo
      if (my_pe==0) print'(/a/)','cancelling MPI service'
      call mpi_finalize(ierr)
end program main


subroutine setup
!=======================================================================
!  setup everything 
!=======================================================================
  use main_module   
  use eke_module   
  use tke_module   
  use idemix_module   
  implicit none
  integer :: i,j,k
  real*8 :: fxa,fxb
 
 if (my_pe==0) print'(/a/)','allocating work space'

!--------------------------------------------------------------
! allocate everything
!--------------------------------------------------------------
  call set_parameter
  call pe_decomposition
  call allocate_main_module
  call allocate_isoneutral_module
  call allocate_tke_module
  call allocate_eke_module
  call allocate_idemix_module
  call allocate_rossmix_module
  call allocate_rossmix2_module
  call allocate_obc_module
  call allocate_tracer_module
  
  call fortran_barrier
  if (my_pe==0) print'(/a/)','done allocating work space'

  if (my_pe==0) print'(/a/)','setting up everything'

!--------------------------------------------------------------
!  Grid
!--------------------------------------------------------------
  call set_grid
  call calc_grid

!--------------------------------------------------------------
! Coriolis
!--------------------------------------------------------------
  call set_coriolis
  call calc_beta

!--------------------------------------------------------------
! topography
!--------------------------------------------------------------
  call set_topography
  call calc_topo
  call calc_spectral_topo

!--------------------------------------------------------------
! other eddy closures
!--------------------------------------------------------------
  call rossmix_initialize
  call rossmix2_initialize  
  call biharmonic_thickness_initialize
!--------------------------------------------------------------
! initial condition and forcing 
!--------------------------------------------------------------
  call set_initial_conditions
  call calc_initial_conditions
  call set_forcing
  if (enable_streamfunction) call streamfunction_init

!--------------------------------------------------------------
! initialize diagnostics
!--------------------------------------------------------------
  call init_diagnostics
  call set_diagnostics
!--------------------------------------------------------------
! initialize EKE module
!--------------------------------------------------------------
  call init_eke
!--------------------------------------------------------------
! initialize isoneutral module
!--------------------------------------------------------------
  call check_isoneutral_slope_crit
!--------------------------------------------------------------
! check setup
!--------------------------------------------------------------
  if (enable_tke .and. .not. enable_implicit_vert_friction) then
    if (my_pe==0) print'(/a)','ERROR: use TKE model only with implicit vertical friction '
    if (my_pe==0) print'(a/)','        -> switch on enable_implicit_vert_fricton        '
    call halt_stop(' in setup')
  endif
  
  if (enable_explicit_vert_friction) then
   fxa=0
   do k=1,nz 
     fxa = max(fxa, kappaM_0*dt_mom*4/dzt(k)**2 )       
   enddo 
   if (my_pe==0) print*,' vertical viscosity max criterium = ',fxa
   if (fxa>1.) then
     if (my_pe==0) print*,'ERROR: kappaM_0 too large, criterium must be <1'
     !call halt_stop('in setup')
   endif 
  endif
  
  if (enable_biharmonic_friction) then
   fxa = 0.
   do j=js_pe,je_pe
    do i=is_pe,ie_pe
     fxb = A_hbi*(cost(j)**biha_friction_cosPower)**2
     fxa = max(fxa, fxb*dt_mom*(4/(cost(j)*dxt(i))**2)**2)
     fxa = max(fxa, fxb*dt_mom*(4/(dyt(j))**2        )**2)
    enddo 
   enddo 
   call global_max(fxa) 
   if (my_pe==0) print*,' biharm. viscosity max criterium = ',fxa
   if (fxa>1.) then
     if (my_pe==0) print*,'ERROR: A_hbi too large, criterium must be <1'
     !call halt_stop('in setup')
   endif
  endif
    
end subroutine setup

