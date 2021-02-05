
!=======================================================================
! Relaxation of a density front
!=======================================================================

module config_module
 !real*8 :: fac=1.0,mix=5e-3
 real*8 :: RES_HOR, RES_TIM
 real*8 :: Ri, alpha
 real*8 :: H0, f0, beta0
 real*8 :: N0, M0, U0 
 real*8 :: Lf, Npcl, Hml
 real*8 :: sigmax, kmax
 real*8 :: dx, dy, dz
 real*8 :: cfl
 real*8 :: betaT = 1.67d-4!, grav=9.81
 real*8 :: lambda_u, lambda_t
 real*8 :: iRe_Av, iRe_Ah, iRe_A4, iPe_Kv, iPe_Kh, iPe_K4
end module config_module



subroutine set_parameter
 use main_module   
 use config_module   
 use diagnostics_module   
 use tke_module   
 use isoneutral_module
 implicit none
 integer :: ntsteps
 character(len=20) :: fmtrea
 character(len=20) :: fmtint
 real*8 :: Ro, delta, L0
 real*8 :: Lx, Ly, Lz 
 real*8 :: para1, para2, para3 
 logical :: fexist

  RES_HOR = 0.5
  RES_TIM = 1

  ! read parameter file
  inquire (file='read_para', exist=fexist)
  if ( fexist ) then
    fmtrea="(ES14.8)"
    open(20, file="read_para", status="old")
    read(20,fmtrea) para1
    read(20,fmtrea) para2
    read(20,fmtrea) para3
    close(20)
  else
    print*, "Warning file read_para does not exist! Take standard values:"
    para1 = 3.16227766E+01  ! Ri
    para2 = 4.0             ! alpha
    para3 = 1.0   
    print*, "para1 = ", para1
    print*, "para2 = ", para2
    print*, "para3 = ", para3
  endif

  !-- start: settings for mixed layer parameterization
  enable_ml_para = .true.
  enable_ml_para_stone       = .true.
  enable_ml_HS_streamfct     = .false.
  ml_algorithm               = 1
  ml_Nsqr_fac                = 4.0
  ! define parameters for mixed layer parameterization
  if ( enable_ml_para_stone) then
    ml_Ce = 0.9
  else 
    ml_Ce = 0.07 ! Fox Kemper 
  endif
  ml_invTau       = 1.6e-6
  !-- end: settings for mixed layer parameterization

  if ( enable_ml_para ) then
    nx  = 1
  else 
    nx  = int(120*RES_HOR)
  endif
  ny    = int(120*RES_HOR)
  nz    = 60
 
  Ri     = para1
  alpha  = para2

  H0     = 300.
  f0     = 7.e-05
  beta0  = 0.0
 
  cfl    = 5.0e-3
 
  N0     = alpha*sqrt(Ri)*f0
  M0     = sqrt(alpha)*f0
  U0     = alpha*H0*f0
 
  sigmax = sqrt(5./54.) / sqrt(1.+Ri) * f0
  kmax   = sqrt(5./2.)  / sqrt(1.+Ri) * f0/U0
 
  dx     = NINT( 2.0*pi/kmax / 10 / RES_HOR) 
  dy     = dx
  dz     = 5.0

  Lf   = ny*dy*0.33 !40.0*dx
  Hml  = 200.
  Npcl = 120*f0 ! large pycnocline stratification
 
  iRe_Ah    = 0.;
  iRe_Av    = 5e-3;
  iRe_A4    = 0.500;
  
  iPe_Kh    = 0;
  iPe_Kv    = iRe_Av;
  iPe_K4    = iRe_A4;
  
  kappaM_0  = U0*dz * iRe_Av
  kappaH_0  = 0.0 !U0*dz * iPe_Kv
  A_hbi     = U0*dx**3 * iRe_A4
  K_hbi     = 0.0 !U0*dx**3 * iPe_K4

  dt_mom    = NINT( cfl * dx/U0 / RES_TIM)
  dt_tracer = dt_mom

  !ntsteps     = 10000
  !runlen      = ntsteps*dt_mom 
  !snapint     = ntsteps*dt_mom/100.
  !ts_monint   = snapint
  ntsteps   = 172800
  runlen    = ntsteps     * RES_HOR*RES_TIM * dt_mom
  ts_monint = ntsteps/800 * RES_HOR*RES_TIM * dt_mom
  !ts_monint = dt_mom
  snapint   = ntsteps/800 * RES_HOR*RES_TIM * dt_mom

  ! parameter for diagnostics
  Ro = 1.0/sqrt(Ri) 
  delta = f0/N0
  L0 = N0*H0/f0
  Lx = nx*dx
  Ly = ny*dy
  Lz = nz*dz

  if (my_pe == 0 ) then
    print*, "================================================================================"
    print*, " Settings for front setup"
    print*, "================================================================================"
    print*, "nx, ny, nz = ", nx, ny, nz
    print*, "Ri = ", Ri, "alpha = ", alpha
    print*, "H0 = ", H0, "f0 = ", f0
    print*, "N0 = ", N0, "M0 = ", M0, "U0 = ", U0
    print*, "Npcl = ", Npcl, "Hml = ", Hml, "Lf = ", Lf
    print*, "dx = ", dx, "dy = ", dy, "dz = ", dz, "dt = ", dt_mom
    print*, "sigmax = ", sigmax, "kmax = ", kmax
    print*, "snapint [days] = ", snapint/86400.0
    print*, "================================================================================"
  endif

  enable_conserve_energy = .false. 
  coord_degree         =.false.
  enable_cyclic_x      =.true.
  enable_hydrostatic   =.true.
  eq_of_state_type = 1 

  congr_epsilon = 1d-12
  congr_max_iterations = 5000
  !enable_streamfunction = .true.

  !congr_epsilon_non_hydro=   1d-6
  !congr_max_itts_non_hydro = 5000    

  enable_explicit_vert_friction = .true. 
  enable_biharmonic_friction    = .true.
  enable_biharmonic_mixing      = .false.
  enable_superbee_advection      = .true.

  enable_diag_ts_monitor = .true.; 
  enable_diag_snapshots  = .true.; 

end subroutine set_parameter


subroutine set_grid
 use main_module   
 use config_module   
 implicit none
 dxt(:) = dx 
 dyt(:) = dy
 dzt(:) = dz 
end subroutine set_grid

subroutine set_coriolis
 use main_module   
 use config_module   
 implicit none
 coriolis_t = f0 
 !coriolis_h = 2*omega*cos( 30./180.*pi)
end subroutine set_coriolis



real*8  function u_star(k)
 use main_module   
 implicit none
 integer :: k
 u_star = 0.6+0.5*tanh( (zt(k)-zt(nz/2))/zt(1)*100)
end function u_star



subroutine set_initial_conditions
  use main_module   
  use config_module
  implicit none
  integer :: i,j,k
  real*8 :: N2(nz)
  real*8 :: dyB(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,1:nz)
  real*8 :: T_perturb(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,1:nz)
  real*8 :: Ly, tv
  logical :: enable_random_perturbation=.true.
  real*8 :: rnum
  !real*8 :: fxa,t_star,u_star
  !real*8 :: betaT = 2.e-4, grav=9.81
  where( zt > -Hml )
    N2 = N0**2
  elsewhere
    N2 = Npcl**2
  end where

  do j=js_pe-onx,je_pe+onx
    do i=is_pe-onx,ie_pe+onx
        if (enable_random_perturbation) then
          call random_number(rnum)
          T_perturb(i,j,:) = 0.001*rnum
        else
          T_perturb(i,j,:) = 0.001*sin(xt(i)*kmax)
        endif
    enddo
  enddo

  Ly = dy*ny
  do j=js_pe-onx,je_pe+onx
    do i=is_pe-onx,ie_pe+onx
      tv = 0.0
      do k=1,nz
        tv = tv + N2(k)*dz
        temp(i,j,k,:) = maskT(i,j,k)*( 1.0/(grav*betaT)*( &
           tv - Lf*M0**2/2.0 * tanh(2.0*(yt(j)-Ly/2.0)/Lf) ) + &
           0.0*T_perturb(i,j,k) + 8.0 )
      enddo
    enddo
  enddo
  u = 0.0
  do j=js_pe-onx,je_pe+onx-1
    do i=is_pe-onx,ie_pe+onx
      do k=1,nz
        dyB(i,j,k) = grav*betaT * (temp(i,j+1,k,tau)-temp(i,j,k,tau))/dyt(j) * maskV(i,j,k)
      enddo
    enddo
  enddo

 do j=js_pe,je_pe
    do i=is_pe,ie_pe
      do k=1,nz-1
        u(i,j,k+1,:) = u(i,j,k,:) - 0.125*( dyB(i,j,k)   + dyB(i+1,j,k)   + dyB(i,j-1,k)   + dyB(i+1,j-1,k)     &
                                          + dyB(i,j,k+1) + dyB(i+1,j,k+1) + dyB(i,j-1,k+1) + dyB(i+1,j-1,k+1) ) &
                                          / coriolis_t(i,j) * dzw(k)
      enddo
    enddo
 enddo

 do j=1,3
    call border_exchg_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,u(:,:,:,j))
    call setcyclic_xyz (is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,u(:,:,:,j))
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
