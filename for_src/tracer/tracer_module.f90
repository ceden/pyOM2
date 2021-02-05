

module tracer_module
     implicit none
!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
     logical :: enable_tracer = .false.
     integer :: n_trac = 0
     real*8, allocatable :: dtrac(:,:,:,:,:)  ! tendency due to advection using Adam Bashforth
     real*8, allocatable :: trac(:,:,:,:,:)   ! tracer
     real*8, allocatable :: trac_source(:,:,:,:)   ! tracer source
     real*8, allocatable :: forc_trac_surface(:,:,:)   ! tracer
end module tracer_module

subroutine allocate_tracer_module
  use main_module
  use tracer_module
  if (enable_tracer) then
   if (n_trac <= 0 ) then
    if (my_pe==0) print*,' ERROR: n_trac <=0 in allocate_tracer_module, set n_trac to positive number larger zero'
    call halt_stop('in allocate_tracer_module')
   endif
   allocate(dtrac(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz,3,n_trac) );dtrac = 0.0
   allocate( trac(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz,3,n_trac) ); trac = 0.0
   allocate( forc_trac_surface(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,n_trac) ); forc_trac_surface = 0.0
   allocate( trac_source(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz,n_trac) ); trac_source = 0.0
  endif
end subroutine allocate_tracer_module


