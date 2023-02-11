




subroutine integrate_tracer
!=======================================================================
! integrate passive tracer
!=======================================================================
 use main_module   
 use isoneutral_module   
 use tracer_module   
 use rossmix_module   
 use timing_module   
 implicit none
 integer :: i,j,k,ks,n
 real*8 :: a_tri(nz),b_tri(nz),c_tri(nz),d_tri(nz),delta(nz)

if (enable_tracer) then
 do n=1,n_trac
  call advect_tracer(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,trac(:,:,:,tau,n),dtrac(:,:,:,tau,n))

  if (enable_AB_time_stepping) then
   !---------------------------------------------------------------------------------
   ! Adam Bashforth time stepping for advection
   !---------------------------------------------------------------------------------
   trac(:,:,:,taup1,n)=trac(:,:,:,tau,n)+dt_tracer*( (1.5+AB_eps)*dtrac(:,:,:,tau,n) - ( 0.5+AB_eps)*dtrac(:,:,:,taum1,n))*maskT
  else
   trac(:,:,:,taup1,n)=trac(:,:,:,tau,n)+dt_tracer*dtrac(:,:,:,tau,n)*maskT
  endif

  !---------------------------------------------------------------------------------
  ! horizontal diffusion
  !---------------------------------------------------------------------------------
  call tic('iso')
  if (enable_hor_diffusion)     call tracer_diffusion(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,trac(:,:,:,:,n) )

  if (enable_biharmonic_mixing) call tracer_biharmonic(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,trac(:,:,:,:,n) )

  !---------------------------------------------------------------------------------
  ! sources like restoring zones, etc
  !---------------------------------------------------------------------------------
  trac(:,:,:,taup1,n)=trac(:,:,:,taup1,n)+dt_tracer*trac_source(:,:,:,n)*maskT

  !---------------------------------------------------------------------------------
  ! isopycnal diffusion
  !---------------------------------------------------------------------------------
  if (enable_neutral_diffusion) then
   call isoneutral_diffusion_trac(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,trac(:,:,:,:,n))
  endif
  if (enable_ml_para) then
   call mixed_layer_para_trac(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,trac(:,:,:,:,n))
  endif

  !---------------------------------------------------------------------------------
  !---------------------------------------------------------------------------------
  !if (enable_rossmix .and. enable_rossmix_bolus_form) then
  if (enable_rossmix) then
    call rossmix_eddy_advect(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,trac(:,:,:,:,n))
  endif

  call toc('iso')

  call tic('vmix')
  !---------------------------------------------------------------------------------
  ! vertical mixing of tracer
  !---------------------------------------------------------------------------------
  a_tri=0.0;b_tri=0.0; c_tri=0.0; d_tri=0.0; delta=0.0
  do j=js_pe,je_pe
   do i=is_pe,ie_pe
    ks=kbot(i,j)
    if (ks>0) then
     do k=ks,nz-1
      delta(k) = dt_tracer/dzw(k)*kappaH(i,j,k)
     enddo
     delta(nz)=0.0
     do k=ks+1,nz
       a_tri(k) = - delta(k-1)/dzt(k)
     enddo
     a_tri(ks)=0.0
     do k=ks+1,nz-1
      b_tri(k) = 1+ delta(k)/dzt(k) + delta(k-1)/dzt(k) 
     enddo
     b_tri(nz) = 1+ delta(nz-1)/dzt(nz) 
     b_tri(ks) = 1+ delta(ks)/dzt(ks)   
     do k=ks,nz-1
      c_tri(k) = - delta(k)/dzt(k)
     enddo
     c_tri(nz)=0.0
     d_tri(ks:nz)=trac(i,j,ks:nz,taup1,n) 
     d_tri(nz) = d_tri(nz) + dt_tracer*forc_trac_surface(i,j,n)/dzt(nz)
     call solve_tridiag(a_tri(ks:nz),b_tri(ks:nz),c_tri(ks:nz),d_tri(ks:nz),trac(i,j,ks:nz,taup1,n),nz-ks+1)
    endif
   enddo
  enddo
  call toc('vmix')

  !---------------------------------------------------------------------------------
  ! open boundaries
  !---------------------------------------------------------------------------------
  call set_obc_tracer(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,trac(:,:,:,:,n),.false.)
 
  !---------------------------------------------------------------------------------
  ! boundary exchange
  !---------------------------------------------------------------------------------
  call border_exchg_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,trac(:,:,:,taup1,n)) 
  call setcyclic_xyz   (is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,trac(:,:,:,taup1,n))
  call set_obc_boundary_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,trac(:,:,:,taup1,n))
  
 enddo
endif
end subroutine integrate_tracer








subroutine write_tracer_restart
!=======================================================================
!=======================================================================
  use main_module   
  use tracer_module   
  implicit none
  character*80 :: filename
  integer :: ierr,io,is,ie,js,je,n
  if (enable_tracer) then
     is=is_pe-onx; ie=ie_pe+onx; js=js_pe-onx; je=je_pe+onx
     write(filename,'(a,i5,a)')  'restart_trac_PE_',my_pe,'.dta'
     call replace_space_zero(filename)
     if (my_pe==0) print'(a,a)',' writing restart tracer file ',filename(1:len_trim(filename))
     call get_free_iounit(io,ierr)
     if (ierr/=0) goto 10
     open(io,file=filename,form='unformatted',status='unknown')
     write(io,err=10) nx,ny,nz
     write(io,err=10) is,ie,js,je
     do n=1,n_trac
       write(io,err=10)  trac(is:ie,js:je,:,tau,n),  trac(is:ie,js:je,:,taum1,n)
       write(io,err=10)  dtrac(is:ie,js:je,:,tau,n),  dtrac(is:ie,js:je,:,taum1,n)
     enddo
     close(io)
     call fortran_barrier()
     return
     10 continue
     print'(a)',' Warning: error writing restart file'
   endif
end subroutine write_tracer_restart



subroutine read_tracer_restart
!=======================================================================
!=======================================================================
  use main_module   
  use tracer_module   
  implicit none
  integer :: ierr,is_,ie_,js_,je_,is,ie,js,je,n
  integer :: io,nx_,ny_,nz_
  character*80 :: filename
  logical :: file_exists
  
  if (enable_tracer) then
     is=is_pe-onx; ie=ie_pe+onx; js=js_pe-onx; je=je_pe+onx
     write(filename,'(a,i5,a)')  'restart_trac_PE_',my_pe,'.dta'
     call replace_space_zero(filename)
     inquire ( FILE=filename, EXIST=file_exists )
     if (.not. file_exists) then
       if (my_pe==0) then
         print'(a,a)',' found no tracer restart file ',filename(1:len_trim(filename))
         print'(a)',' proceeding with initial conditions'
       endif
       return
     endif

     if (my_pe==0) print'(a,a)',' reading tracer from restart file ',filename(1:len_trim(filename))
     call get_free_iounit(io,ierr)
     if (ierr/=0) goto 10
     open(io,file=filename,form='unformatted',status='old',err=10)
     read(io,err=10) nx_,ny_,nz_
     if (nx/=nx_ .or. ny/=ny_ .or. nz/= nz_) then 
       if (my_pe==0) then
        print*,' read from restart dimensions: ',nx_,ny_,nz_
        print*,' does not match dimensions   : ',nx,ny,nz
       endif
       goto 10
     endif
     read(io,err=10) is_,ie_,js_,je_
     if (is_/=is.or.ie_/=ie.or.js_/=js.or.je_/=je) then
       if (my_pe==0) then
        print*,' read from restart PE boundaries: ',is_,ie_,js_,je_
        print*,' which does not match           : ',is,ie,js,je
       endif
       goto 10
     endif

     do n=1,n_trac
       read(io,err=10)  trac(is:ie,js:je,:,tau,n), trac(is:ie,js:je,:,taum1,n)
       read(io,err=10) dtrac(is:ie,js:je,:,tau,n),dtrac(is:ie,js:je,:,taum1,n)
     enddo
     close(io)
     call fortran_barrier()
   
     return
     10 continue
     print'(a)',' Warning: error reading tracer restart file'
     call halt_stop(' in read_restart')
   endif
end subroutine read_tracer_restart

