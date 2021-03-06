
!=======================================================================
!  test for isopycnal mixing
!======================================================================= 


module config_module
 ! use this module only locally in this file
 implicit none
 real*8 :: y0 = -40, y1 = 28-60, y2=32-60
 real*8 :: hresolv = 5.0, vresolv = 2.0, yt_1,yt_ny,xt_1,xt_nx
end module config_module
 
 

subroutine set_parameter
 ! ----------------------------------
 !       set here main parameter
 ! ----------------------------------
 use main_module   
 use eke_module   
 use tke_module   
 use idemix_module   
 use isoneutral_module   
 use diagnostics_module   
 use config_module
 implicit none

  nx   = int(20*hresolv); nz   = int(10*vresolv); ny  = int(20*hresolv)
  dt_mom    = 1200.0
  dt_tracer = 1200.0

  coord_degree     = .true.
  enable_cyclic_x  = .true.

  runlen = 365*86400.*1  

 enable_diag_ts_monitor = .true.; ts_monint = dt_tracer
 enable_diag_snapshots  = .true.; snapint  =  86400.
 enable_diag_energy  = .true.; energfreq =  dt_tracer; energint  =  dt_tracer*10
 enable_diag_variances = .true.; varfreq =  dt_tracer; varint  =  dt_tracer*100

 enable_diag_averages   = .true.
 aveint  = 365*86400
 avefreq = dt_tracer*10

  congr_epsilon = 1e-9
  enable_streamfunction = .true.
  congr_max_iterations = 15000

  enable_bottom_friction = .true.; r_bot = 1e-5*vresolv

  enable_conserve_energy = .true.
  enable_store_bottom_friction_tke = .true.
  enable_store_cabbeling_heat = .true.
  enable_take_P_diss_adv_from_tke = .false.

 !enable_biharmonic_mixing = .true.
 !K_hbi  = 0e11!/hresolv**4
 enable_superbee_advection = .true.

  enable_biharmonic_friction  = .true. ! for eddy resolving version
  A_hbi  = 5e11!/hresolv**4

  enable_implicit_vert_friction = .true.; 
  enable_tke = .true.
  c_k = 0.1
  c_eps = 0.7
  alpha_tke = 30.0
  mxl_min = 1d-8
  tke_mxl_choice = 2
  enable_tke_superbee_advection = .false.

  enable_idemix = .true.
  enable_idemix_hor_diffusion = .true.; 
  enable_idemix_superbee_advection = .true.

  eq_of_state_type = 5 
end subroutine set_parameter



subroutine set_grid
  use main_module   
  use config_module   
  implicit none
  dxt = 1.0/hresolv
  dyt = 1.0/hresolv
  x_origin= 0.0
  y_origin= y0
  dzt = 1000.0/nz

 yt_1  = y0-dyt(js_pe)/2.0
 yt_ny = yt_1 + dyt(js_pe)*ny
 xt_1  = 0.-dxt(is_pe)/2.0
 xt_nx = xt_1 + dxt(is_pe)*nx
end subroutine set_grid


subroutine set_coriolis
 use main_module   
 use config_module   
 implicit none
 real*8 :: phi0,betaloc
 integer :: j
 phi0 = (y2+y1)/2.0 /180. *pi
 betaloc = 2*omega*cos(phi0)/radius
 do j=js_pe-onx,je_pe+onx
  coriolis_t(:,j) = 2*omega*sin(phi0) +betaloc*yt(j)
 enddo

end subroutine set_coriolis


subroutine set_initial_conditions
 use main_module   
 use config_module   
 use idemix_module   
 use tke_module   
 implicit none
 integer :: i,j,k
 real*8 :: y
 do k=1,nz
  do j=js_pe-onx,je_pe+onx
   do i=is_pe-onx,ie_pe+onx
     y = yt(j)/yt_ny
    salt(i,j,k,:) = 35
    temp(i,j,k,:) =  (1-zt(k)/zw(1))*15 + &
        2*(tanh((y-0.5)/0.1)+1.0) +  3*(tanh(-(y-0.5)/0.2)+1.0)
   enddo
  enddo
 enddo


 do j=js_pe-onx,je_pe+onx
   do i=is_pe-onx,ie_pe+onx
     if ( yt(j)< y1) surface_taux(i,j) = &
         .1e-3*sin(pi*(yu(j)-y0)/(y1-yt_1))*maskU(i,j,nz)
   enddo
 enddo

 if (enable_idemix ) then
  do j=js_pe-onx,je_pe+onx
   do i=is_pe-onx,ie_pe+onx
     forc_iw_bottom(i,j) =  3e-6
     forc_iw_surface(i,j) = 0.5e-6
   enddo
  enddo
 endif

 if (enable_tke ) then
  do j=js_pe-onx+1,je_pe+onx
   do i=is_pe-onx+1,ie_pe+onx
     forc_tke_surface(i,j) = sqrt( (0.5*(surface_taux(i,j)+surface_taux(i-1,j)))**2  &
                                  +(0.5*(surface_tauy(i,j)+surface_tauy(i,j-1)))**2 )**(3./2.) 
   enddo
  enddo
 endif

end subroutine set_initial_conditions




function tstar(j)
 use main_module   
 use config_module   
 implicit none
 integer :: j
 real*8 :: tstar
 tstar=15
 if (yt(j)<y1) tstar=15*(yt(j)-yt_1)/(y1-yt_1)
 if (yt(j)>y2) tstar=15*(1-(yt(j)-y2)/(yt_ny-y2) )
end function tstar



subroutine set_forcing
 use main_module   
 implicit none
 integer :: i,j
 real*8 :: tstar

 do j=js_pe-onx,je_pe+onx
  do i=is_pe-onx,ie_pe+onx
    forc_temp_surface(i,j)=dzt(nz)/(30.*86400.)*(tstar(j)-temp(i,j,nz,tau)) 
  enddo
 enddo
end subroutine set_forcing


subroutine set_topography
 use main_module   
 use config_module   
 implicit none
 integer :: i,j
 !kbot=0
 kbot=1
 do i=is_pe,ie_pe
   do j=js_pe,je_pe
     !if ( yt(j)>=y1.and.xt(i)<=int(xt_nx*0.1) )  kbot(i,j)=0
     if ( yt(j)>=y1.and.i>=nx/2-hresolv.and.i<=nx/2+hresolv)  kbot(i,j)=0
   enddo
 enddo
end subroutine set_topography



subroutine set_diagnostics


 use main_module   
 use eke_module   
 use tke_module   
 use idemix_module   
 use isoneutral_module   
 implicit none
 call register_average('taux','Zonal wind stress','m^2/s','UT',surface_taux,0D0,.false.)
 call register_average('tauy','Meridional wind stress','m^2/s','TU',surface_tauy,0D0,.false.)
 call register_average('forc_temp_surface','Surface temperature flux','m K/s','TT',forc_temp_surface,0D0,.false.)
 call register_average('forc_salt_surface','Surface salinity flux','m g/s kg','TT',forc_salt_surface,0D0,.false.)
 if (enable_streamfunction) then
   call register_average('psi','Barotropic streamfunction','m^2/s','UU',psi(:,:,tau),0D0,.false.)
 else
   call register_average('psi','Surface pressure','m^2/s','TT',psi(:,:,tau),0D0,.false.)
 endif
 call register_average('temp','Temperature','deg C','TTT',0d0,temp(:,:,:,tau),.true.)
 call register_average('salt','Salinity','g/kg','TTT',0d0,salt(:,:,:,tau),.true.)
 call register_average('u','Zonal velocity','m/s','UTT',0d0,u(:,:,:,tau),.true.)
 call register_average('v','Meridional velocity','m/s','TUT',0d0,v(:,:,:,tau),.true.)
 call register_average('w','Vertical velocity','m/s','TTU',0d0,w(:,:,:,tau),.true.)
 call register_average('Nsqr','Square of stability frequency','1/s^2','TTU',0d0,Nsqr(:,:,:,tau),.true.)
 call register_average('Hd','Dynamic enthalpy','m^2/s^2','TTT',0d0,Hd(:,:,:,tau),.true.)
 call register_average('rho','Density','kg/m^3','TTT',0d0,rho(:,:,:,tau),.true.)

 call register_average('K_diss_v','Dissipation by vertical friction','m^2/s^3','TTU',0d0,K_diss_v,.true.)
 call register_average('K_diss_h','Dissipation by lateral friction','m^2/s^3','TTU',0d0,K_diss_h,.true.)
 call register_average('K_diss_bot','Dissipation by bottom friction','m^2/s^3','TTU',0d0,K_diss_bot,.true.)
 call register_average('P_diss_v','Dissipation by vertical mixing','m^2/s^3','TTU',0d0,P_diss_v,.true.)
 call register_average('P_diss_nonlin','Dissipation by nonlinear vert. mix.','m^2/s^3','TTU',0d0,P_diss_nonlin,.true.)

 call register_average('kappaH','Vertical diffusivity','m^2/s','TTU',0d0,kappaH,.true.)

 if (enable_tke)  then
   call register_average('TKE','Turbulent kinetic energy','m^2/s^2','TTU',0d0,tke(:,:,:,tau),.true.)
   call register_average('Prandtl','Prandtl number',' ','TTU',0d0,Prandtlnumber,.true.)
   call register_average('mxl','Mixing length',' ','TTU',0d0,mxl,.true.)
   call register_average('tke_diss','Dissipation of TKE','m^2/s^3','TTU',0d0,tke_diss,.true.)
   call register_average('forc_tke_surface','TKE surface forcing','m^3/s^2','TT',forc_tke_surface,0D0,.false.)
   call register_average('tke_surface_corr','TKE surface flux correction','m^3/s^2','TT',tke_surf_corr,0D0,.false.)
 endif

 if (enable_idemix)  then
   call register_average('E_iw','Internal wave energy','m^2/s^2','TTU',0d0,e_iw(:,:,:,tau),.true.)
   call register_average('forc_iw_surface','IW surface forcing','m^3/s^2','TT',forc_iw_surface,0D0,.false.)
   call register_average('forc_iw_bottom','IW bottom forcing','m^3/s^2','TT',forc_iw_bottom,0D0,.false.)
   call register_average('iw_diss','Dissipation of E_iw','m^2/s^3','TTU',0d0,iw_diss,.true.)
   call register_average('c0','Vertical IW group velocity','m/s','TTU',0d0,c0,.true.)
   call register_average('v0','Horizontal IW group velocity','m/s','TTU',0d0,v0,.true.)
 endif



end subroutine set_diagnostics

subroutine set_particles
end subroutine set_particles
