



subroutine init_diagnostics
!=======================================================================
! initialize diagnostic routines
!=======================================================================
 use main_module
 use tke_module
 use eke_module
 use idemix_module
 use rossmix_module
 use isoneutral_module
 use diagnostics_module   
 use tracer_module   
 implicit none

 if (my_pe==0) print'(/,a)','Diagnostic setup:'

 if (enable_diag_ts_monitor) then
   if (my_pe==0) print'(a,e12.6,a,f10.2,a)',' time step monitor every ',ts_monint,' seconds/',ts_monint/dt_tracer,' time steps'
 endif

 if (enable_diag_tracer_content) then
   if (my_pe==0) print'(a,e12.6,a,f10.2,a)',' monitor tracer content every ',&
                  trac_cont_int,' seconds/',trac_cont_int/dt_tracer,' time steps'
 endif

 if (enable_diag_snapshots) then
    if (my_pe==0) print'(a,e12.6,a,f10.2,a)',' writing snapshots every ',snapint,' seconds/',snapint/dt_tracer,' time steps'
    call init_snap_cdf
    if (enable_tke) call init_snap_cdf_tke
    if (enable_eke) call init_snap_cdf_eke
    if ( (enable_neutral_diffusion .and. enable_skew_diffusion) .or. enable_TEM_friction &
                 .or. enable_ml_para) call init_snap_cdf_isoneutral
    if (enable_conserve_energy) call init_snap_cdf_energy
    if (enable_idemix) call init_snap_cdf_idemix
    if (enable_rossmix) call init_snap_cdf_rossmix
    if (enable_tracer) call init_snap_cdf_tracer
 endif

 if ( enable_diag_parallel_snap )  call init_diag_snap_parallel
 if ( enable_diag_chunks_snap )    call init_diag_snap_chunks
 
 if (enable_diag_averages)  then
    if (my_pe==0) print'(a,e12.6,a,f10.2,a)',' writing time averages every ',aveint,' seconds/',aveint/dt_tracer,' time steps'
    if (my_pe==0) print'(a,f10.2,a)',' calculating every ',avefreq/dt_tracer,' time step'
 endif

 if (enable_diag_variances)  then
    if (my_pe==0) print'(a,e12.6,a,f10.2,a)',' writing variances every ',varint,' seconds/',varint/dt_tracer,' time steps'
    if (my_pe==0) print'(a,f10.2,a)',' calculating every ',varfreq/dt_tracer,' time step'
    call init_diag_variance
 endif

 if (enable_diag_energy) then
    if (my_pe==0) print'(a,e12.6,a,f10.2,a)',' writing energetics every ',energint,' seconds/',energint/dt_tracer,' time steps'
    if (my_pe==0) print'(a,f10.2,a)',' calculating every ',energfreq/dt_tracer,' time step'
    call init_diag_energy
 endif

 if (enable_diag_overturning) then
    if (my_pe==0) print'(a,e12.6,a,f10.2,a)',' writing isopyc. overturning every ', &
                overint,' seconds/',overint/dt_tracer,' time steps'
    if (my_pe==0) print'(a,f10.2,a)',' calculating every ',overfreq/dt_tracer,' time step'
    call init_diag_overturning
 endif

 if (enable_diag_particles) then
    if (my_pe==0) print'(a,e12.6,a,f10.2,a)',' writing particles every ', &
        particles_int,' seconds/',particles_int/dt_tracer,' time steps'
    call set_particles()
    call init_diag_particles
    call init_write_particles
 endif

 if (enable_diag_snap_fluxes) then
   if (my_pe==0) print'(a,e12.6,a,f10.2,a)',' writing fluxes every ', &
        snapint,' seconds/',snapint/dt_tracer,' time steps'
   call init_diag_snap_fluxes
 endif
 
 if (enable_diag_snap_tendency) then
   if (my_pe==0) print'(a,e12.6,a,f10.2,a)',' writing tendencies every ', &
        snapint,' seconds/',snapint/dt_tracer,' time steps'
   call init_diag_snap_tendency
 endif
 
end subroutine init_diagnostics




subroutine diagnose
!=======================================================================
! call diagnostic routines
!=======================================================================
 use main_module   
 use diagnostics_module   
 use tke_module   
 use eke_module   
 use idemix_module   
 use rossmix_module   
 use isoneutral_module   
 use tracer_module   
 implicit none
 logical :: GM_strfct_diagnosed
 real*8 :: time

 GM_strfct_diagnosed = .false.
 time = itt*dt_tracer


 !if ( enable_diag_ts_monitor .and.  modulo(time,ts_monint) < dt_tracer ) then
 if ( enable_diag_ts_monitor .and.  (mod(itt,int(ts_monint/dt_tracer)) ==0 .or. itt==0) ) then
    if (my_pe==0 .and. enable_hydrostatic) print'(a,i10.10,a,e10.4,a,i6)', &
                      ' itt=',itt,' time=',time,'s solver itts=',congr_itts
    if (my_pe==0 .and. .not. enable_hydrostatic) print'(a,i10.10,a,e10.4,a,i6,a,i6,a)', &
                      ' itt=',itt,' time=',time,'s solver itts=',congr_itts,'( non hydro itts= ',congr_itts_non_hydro,')'
   call diag_cfl
 endif


 !if ( enable_diag_snapshots.and.  modulo(time,snapint) < dt_tracer  )   then
 if ( enable_diag_snapshots .and.  (mod(itt,int(snapint/dt_tracer)) ==0 .or. itt==0) ) then
    if (enable_neutral_diffusion .and. enable_skew_diffusion .and. .not.GM_strfct_diagnosed) then
      call isoneutral_diag_streamfunction
      GM_strfct_diagnosed = .true.
    endif
    call diag_snap
    if (enable_tke) call diag_snap_tke
    if (enable_eke) call diag_snap_eke
    if ( (enable_neutral_diffusion .and. enable_skew_diffusion) .or. enable_TEM_friction &
               .or. enable_ml_para) call diag_snap_isoneutral
    if (enable_conserve_energy) call diag_snap_energy
    if (enable_idemix) call diag_snap_idemix
    if (enable_rossmix) call diag_snap_rossmix
    if (enable_tracer) call diag_snap_tracer
 endif

 if ( enable_diag_parallel_snap .and.  (mod(itt,int(snapint/dt_tracer)) ==0 .or. itt==0) ) then
   call diag_snap_parallel
 endif 

 if ( enable_diag_chunks_snap .and.  (mod(itt,int(snapint/dt_tracer)) ==0 .or. itt==0) ) then
   call diag_snap_chunks
 endif 

 if ( enable_diag_tracer_content .and.  modulo(time,trac_cont_int) < dt_tracer ) then
   call diag_tracer_content
 endif


 if ( enable_diag_energy .and.  modulo(time,energfreq) < dt_tracer  )  call diagnose_energy
 if ( enable_diag_energy .and.  modulo(time,energint)  < dt_tracer  )  call write_energy

 if ( enable_diag_averages .and.  modulo(time,avefreq) < dt_tracer  ) then
    if (enable_neutral_diffusion .and. enable_skew_diffusion .and. .not.GM_strfct_diagnosed) then
      call isoneutral_diag_streamfunction
      GM_strfct_diagnosed = .true.
    endif
    call diag_averages
 endif

 if ( enable_diag_averages .and.  modulo(time,aveint) < dt_tracer .and. itt>1 )  call write_averages

 if ( enable_diag_variances .and.  modulo(time,varfreq) < dt_tracer  ) call diag_variance
 if ( enable_diag_variances .and.  modulo(time,varint)  < dt_tracer .and. itt>1 )  call diag_variance_write


 if ( enable_diag_overturning .and.  modulo(time,overfreq) < dt_tracer )  then
    if (enable_neutral_diffusion .and. enable_skew_diffusion .and. .not.GM_strfct_diagnosed) then
      call isoneutral_diag_streamfunction
      GM_strfct_diagnosed = .true.
    endif
    call diag_overturning
 endif

 if ( enable_diag_overturning .and.  modulo(time,overint) < dt_tracer )  call write_overturning


 if ( enable_diag_particles)   then
    call integrate_particles
    if ( modulo(time,particles_int) < dt_tracer  )   call write_particles
 endif
 
 
 if ( enable_diag_snap_fluxes .and.  (mod(itt,int(snapint/dt_tracer)) ==0 .or. itt==0) ) then
   call diag_snap_fluxes
 endif 

 if ( enable_diag_snap_tendency .and.  (mod(itt,int(snapint/dt_tracer)) ==0 .or. itt==0) ) then
   call diag_snap_tendency
 endif 

end subroutine diagnose




subroutine diag_cfl
!=======================================================================
! check for CFL violation
!=======================================================================
 use main_module   
 use tke_module   
 use eke_module   
 use idemix_module   
 use diagnostics_module   
 implicit none
 integer :: i,j,k,ii,jj,kk
 real*8 :: cfl,wcfl

 cfl = 0d0; wcfl=0d0
 do k=1,nz
  do j=js_pe,je_pe
   do i=is_pe,ie_pe
     cfl = max(cfl,  abs(u(i,j,k,tau))*maskU(i,j,k)/(cost(j)*dxt(i))*dt_tracer )
     cfl = max(cfl,  abs(v(i,j,k,tau))*maskV(i,j,k)/dyt(j)*dt_tracer )
     wcfl = max(wcfl,  abs(w(i,j,k,tau))*maskW(i,j,k)/dzt(k)*dt_tracer )
   enddo
  enddo
 enddo
 call global_max(cfl); call global_max(wcfl)
 
 if  ((     enable_hydrostatic .and. (cfl > 1.0) ) .or. &
      (.not.enable_hydrostatic .and. (cfl > 1.0 .or. wcfl > 1.0) )) then
   if (my_pe==0) print'(/a)','ERROR:  CFL number violation'
   if (my_pe==0) print'(a,f12.6)','max. horizontal CFL number = ',cfl
   if (my_pe==0) print'(a,f12.6)','max. vertical CFL number   = ',wcfl
   ii=-1;jj=-1;kk=-1
   do k=1,nz
    do j=js_pe,je_pe
     do i=is_pe,ie_pe
      if (cfl <= abs(u(i,j,k,tau))*maskU(i,j,k)/(cost(j)*dxt(i))*dt_tracer  .or. &
          cfl <=  abs(v(i,j,k,tau))*maskV(i,j,k)/dyt(j)*dt_tracer     .or. &
          (wcfl <= abs(w(i,j,k,tau))*maskW(i,j,k)/dzt(k)*dt_tracer .and. .not.enable_hydrostatic)  ) then
        ii=i;jj=j;kk=k
      endif
     enddo
    enddo
   enddo 
   call global_max_int(ii)
   call global_max_int(jj)
   call global_max_int(kk)
   if (my_pe==0) print'(a,i5,a,i5,a,i5)', 'max CFL number at i=',ii,' j=',jj,' k=',kk
   if (my_pe==0) print'(a,i9,a/)' ,' at itt = ',itt,' ... stopping integration '
   call panic_snap
   call halt_stop(' in diag_cfl')
 endif
 
 ! check for NaN
 if (cfl/=cfl .or. wcfl /=wcfl) then
   if (my_pe==0) print'(/a)','ERROR:   CFL number is NaN '
   if (my_pe==0) print'(a,i9,a/)' ,' at itt = ',itt,' ... stopping integration '
   call panic_snap
   call halt_stop(' in diag_cfl')
 endif
 
 if (my_pe==0) print'(a,f12.6,f12.6)', ' maximal hor./ver CFL number =',cfl,wcfl

 !if (enable_eke .or. enable_tke .or. enable_idemix) then
 ! cfl = 0d0; wcfl=0d0
 ! do k=1,nz
 !  do j=js_pe,je_pe
 !   do i=is_pe,ie_pe
 !     cfl = max(cfl,  abs(u_wgrid(i,j,k))*maskU(i,j,k)/(cost(j)*dxt(i))*dt_tracer )
 !     cfl = max(cfl,  abs(v_wgrid(i,j,k))*maskV(i,j,k)/dyt(j)*dt_tracer )
 !     wcfl = max(wcfl,  abs(w_wgrid(i,j,k))*maskW(i,j,k)/dzt(k)*dt_tracer )
 !   enddo
 !  enddo
 ! enddo
 ! call global_max(cfl); call global_max(wcfl)
 ! !if (my_pe==0) print'(a,f12.6)', ' maximal hor. CFL number on w grid=',cfl
 ! !if (my_pe==0) print'(a,f12.6)', ' maximal ver. CFL number on w grid=',wcfl
 !endif

end subroutine diag_cfl






subroutine diag_tracer_content
!=======================================================================
! Diagnose tracer content
!=======================================================================
 use main_module   
 implicit none
 integer :: i,j,k
 real*8 :: fxa,tempm,saltm,volm,vtemp,vsalt
 real*8, save :: tempm1=0.,saltm1=0.,vtemp1=0.,vsalt1=0.

  volm=0;tempm=0;vtemp=0;saltm=0;vsalt=0
  do k=1,nz 
   do j=js_pe,je_pe
    do i=is_pe,ie_pe
     fxa = area_t(i,j)*dzt(k)*maskT(i,j,k)
     volm = volm + fxa
     tempm = tempm + fxa*temp(i,j,k,tau)
     saltm = saltm + fxa*salt(i,j,k,tau)
     vtemp = vtemp + temp(i,j,k,tau)**2*fxa
     vsalt = vsalt + salt(i,j,k,tau)**2*fxa
    enddo
   enddo
  enddo
  call global_sum(tempm); call global_sum(saltm); call global_sum(volm);  
  call global_sum(vtemp); call global_sum(vsalt);   

  if (my_pe==0) then
     print*,'  '
     print'(a,f20.15,a,e20.14)',' mean temperature ',tempm /volm,' change to last ',(tempm-tempm1)/volm
     print'(a,f20.15,a,e20.14)',' mean salinity    ',saltm /volm,' change to last ',(saltm-saltm1)/volm
     print'(a,e20.14,a,e20.14)',' temperature var. ',vtemp /volm,' change to last ',(vtemp-vtemp1)/volm
     print'(a,e20.14,a,e20.14)',' salinity var.    ',vsalt /volm,' change to last ',(vsalt-vsalt1)/volm
  endif

  tempm1=tempm; vtemp1=vtemp; saltm1=saltm; vsalt1=vsalt
end subroutine diag_tracer_content




subroutine panic_snap
!=======================================================================
! write a complete snapshot before panic shutdown
!=======================================================================
 use main_module
 use tke_module
 use eke_module
 use idemix_module
 use rossmix_module
 use isoneutral_module
 use diagnostics_module
 implicit none
 
 if (.not. enable_diag_snapshots) then
    call init_snap_cdf
    if (enable_tke) call init_snap_cdf_tke
    if (enable_eke) call init_snap_cdf_eke
    if ( (enable_neutral_diffusion .and. enable_skew_diffusion) .or. enable_TEM_friction) call init_snap_cdf_isoneutral
    if (enable_conserve_energy) call init_snap_cdf_energy
    if (enable_idemix) call init_snap_cdf_idemix
    if (enable_rossmix) call init_snap_cdf_rossmix
 endif

 if (enable_neutral_diffusion .and. enable_skew_diffusion) call isoneutral_diag_streamfunction
 call diag_snap
 if (enable_tke) call diag_snap_tke
 if (enable_eke) call diag_snap_eke
 if ( (enable_neutral_diffusion .and. enable_skew_diffusion) .or. enable_TEM_friction) call diag_snap_isoneutral
 if (enable_conserve_energy) call diag_snap_energy
 if (enable_idemix) call diag_snap_idemix
 if (enable_rossmix) call diag_snap_rossmix
 
end subroutine panic_snap

