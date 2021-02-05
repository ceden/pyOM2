






program time_step
!=======================================================================
!     just do one time step for Machenhauer method
!=======================================================================
      use main_module   
      implicit none
      integer :: ierr,iargc
      character (len=80) :: arg
!---------------------------------------------------------------------------------
!       Initialize MPI setting
!---------------------------------------------------------------------------------
      call mpi_init(ierr)
      call my_mpi_init(my_comm)

      n_pes_i = 1; n_pes_j = 1
      if (n_pes>1) then
       if (iargc() < 2) then
        call halt_stop(' not enough command line input')
       endif
       call getarg(1,arg); read(arg,*) n_pes_i
       call getarg(2,arg); read(arg,*) n_pes_j
       if (my_pe==0) print'(/a,i4,a,i4,a)','using ',n_pes_i,' x ',n_pes_j ,' PEs'
      endif
!--------------------------------------------------------------
!     allocate everything
!--------------------------------------------------------------
      call set_parameter
      call pe_decomposition
      call allocate_main_module
!--------------------------------------------------------------
!      Grid
!--------------------------------------------------------------
      call set_grid
      call calc_grid
!--------------------------------------------------------------
!     topography
!--------------------------------------------------------------
      kbot=1
      call calc_topo
!--------------------------------------------------------------
!     initial condition and forcing 
!--------------------------------------------------------------
      call set_initial_conditions
      call calc_initial_conditions

      call tendencies


      call mpi_finalize(ierr)
end program time_step






subroutine tendencies
 use main_module
 implicit none
 include "netcdf.inc"
 integer :: ncid,iret,id,xdim,ydim,zdim,npe
 real*8 :: alpha,get_drhodT


 call momentum_advection
 call advect_tracer(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,temp(:,:,:,tau),dtemp(:,:,:,tau))
 alpha = get_drhodT(35d0,5d0,0d0) 
 ! T=-B/grav/alpha*rho_0
 dtemp = -dtemp*grav*alpha/rho_0

 if (my_pe==0)  then
  ncid = nccre ('machenhauer1.cdf', NCCLOB, iret)
  iret=nf_set_fill(ncid, NF_NOFILL, iret)
  xdim  = ncddef(ncid, 'x', nx, iret)
  ydim  = ncddef(ncid, 'y', ny, iret)
  zdim  = ncddef(ncid, 'z', nz, iret)
  id  = ncvdef (ncid,'du', NCFLOAT,3,(/xdim,ydim,zdim/),iret)
  id  = ncvdef (ncid,'dv', NCFLOAT,3,(/xdim,ydim,zdim/),iret)
  id  = ncvdef (ncid,'db', NCFLOAT,3,(/xdim,ydim,zdim/),iret)
  call ncclos (ncid, iret)
 endif

 do npe = 0,n_pes-1
  if (my_pe==npe) then
   iret=nf_open('machenhauer1.cdf',NF_WRITE,ncid)
   iret=nf_inq_varid(ncid,'du',id)
   iret= nf_put_vara_double(ncid,id,(/is_pe,js_pe,1/),(/i_blk,j_blk,nz/),du_adv(is_pe:ie_pe,js_pe:je_pe,:))
   iret=nf_inq_varid(ncid,'dv',id)
   iret= nf_put_vara_double(ncid,id,(/is_pe,js_pe,1/),(/i_blk,j_blk,nz/),dv_adv(is_pe:ie_pe,js_pe:je_pe,:))
   iret=nf_inq_varid(ncid,'db',id)
   iret= nf_put_vara_double(ncid,id,(/is_pe,js_pe,1/),(/i_blk,j_blk,nz/),dtemp(is_pe:ie_pe,js_pe:je_pe,:,tau))
   call ncclos (ncid, iret)
  endif
  call fortran_barrier()
 enddo

end subroutine tendencies

