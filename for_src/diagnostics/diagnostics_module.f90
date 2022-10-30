

module diagnostics_module
      implicit none
!---------------------------------------------------------------------------------
!     diagnostic options
!---------------------------------------------------------------------------------
      logical :: enable_diag_ts_monitor    = .false. ! enable time step monitor
      logical :: enable_diag_energy        = .false. ! enable diagnostics for energy
      logical :: enable_diag_averages      = .false. ! enable time averages
      logical :: enable_diag_variances     = .false. ! enable time variances
      logical :: enable_diag_snapshots     = .false. ! enable snapshots
      logical :: enable_diag_overturning   = .false. ! enable isopycnal overturning diagnostic
      logical :: enable_diag_tracer_content= .false. ! enable tracer content and variance monitor
      logical :: enable_diag_particles     = .false. ! enable integration of particles
      logical :: enable_diag_qg_streamfct  = .false. ! enable quasi-geostrophic streamfct. diagnostics
      logical :: enable_diag_snap_fluxes   = .false. ! enable advective and diffusive flux diagnostics
      logical :: enable_diag_snap_tendency = .false. ! enable temp/salt tendency diagnostics
      logical :: enable_diag_parallel_snap = .false. ! enable parallel output
      logical :: enable_diag_chunks_snap   = .false. ! enable snapshot output in chunks per processor
      
      real*8  :: snapint=0.  ! intervall between snapshots to be written in seconds
      real*8  :: aveint=0.   ! intervall between time averages to be written in seconds
      real*8  :: varint=0.   ! intervall between time variances to be written in seconds
      real*8  :: energint=0. ! intervall between energy diag to be written in seconds
      real*8  :: energfreq=0.! diagnosing every energfreq seconds 
      real*8  :: ts_monint=0.! intervall between time step monitor in seconds
      real*8  :: avefreq=0.  ! averaging every ave_freq seconds 
      real*8  :: varfreq=0.  ! calculating variances every var_freq seconds 
      real*8  :: overint=0.  ! intervall between overturning averages to be written in seconds
      real*8  :: overfreq=0. ! averaging overturning every ave_freq seconds 
      real*8  :: trac_cont_int=0.! intervall between tracer content monitor in seconds
      real*8  :: particles_int=0. ! intervall  
end module diagnostics_module





