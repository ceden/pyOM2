


module diag_variance_module
!=======================================================================
! Module for variances
!=======================================================================
  implicit none
  integer :: nitts=0

      real*8, allocatable :: umean(:,:,:),vmean(:,:,:),rhom(:,:,:)
      real*8, allocatable :: uvar(:,:,:),vvar(:,:,:),wmean(:,:,:)
      real*8, allocatable :: urho(:,:,:),vrho(:,:,:),wrho(:,:,:)
      real*8,parameter :: epsln     = 1.0d-20
      real*8,parameter :: spval = -1.0d33
  
end module diag_variance_module




subroutine init_diag_variance
!=======================================================================
!=======================================================================
      use main_module   
      use diag_variance_module
      implicit none

      if (my_pe==0) print*,' Initializing diagnosis of variances'
      nitts=0
      
      allocate(umean(is_pe-onx:ie_pe+onx ,js_pe-onx:je_pe+onx,nz))
      allocate(vmean(is_pe-onx:ie_pe+onx ,js_pe-onx:je_pe+onx,nz))
      allocate(wmean(is_pe-onx:ie_pe+onx ,js_pe-onx:je_pe+onx,nz))
      umean=0; vmean=0; wmean=0
      allocate(uvar(is_pe-onx:ie_pe+onx ,js_pe-onx:je_pe+onx,nz))
      allocate(vvar(is_pe-onx:ie_pe+onx ,js_pe-onx:je_pe+onx,nz))
      uvar=0; vvar=0; 
      allocate(urho(is_pe-onx:ie_pe+onx ,js_pe-onx:je_pe+onx,nz))
      allocate(vrho(is_pe-onx:ie_pe+onx ,js_pe-onx:je_pe+onx,nz))
      allocate(wrho(is_pe-onx:ie_pe+onx ,js_pe-onx:je_pe+onx,nz))
      allocate(rhom(is_pe-onx:ie_pe+onx ,js_pe-onx:je_pe+onx,nz))
      urho=0; vrho=0; wrho=0; rhom=0
      if (my_pe==0) print*,' done initializing diagnosis of variances'
end subroutine init_diag_variance





subroutine diag_variance
!=======================================================================
  !=======================================================================
  ! bar( u (r_i+1 + r_i)  ) = bar( (um + u') (rm_+ + r'_+ + rm_i + r'_i) )
  !                         = um(rm_+ + rm_i) 
  !                          +bar(  u' (r'_+ r'_i) )
      use main_module   
      use diag_variance_module
      implicit none
      integer :: i,j,k
      nitts = nitts + 1
      umean = umean + u(:,:,:,tau)
      vmean = vmean + v(:,:,:,tau)
      wmean = wmean + w(:,:,:,tau)
      uvar  = uvar + u(:,:,:,tau)**2
      vvar  = vvar + v(:,:,:,tau)**2
      rhom  = rhom + rho(:,:,:,tau)
         
      do i=is_pe-onx,ie_pe+onx-1
       urho(i,:,:) = urho(i,:,:)+u(i,:,:,tau)*0.5*( rho(i,:,:,tau)+rho(i+1,:,:,tau) )
      enddo
      do j=js_pe-onx,je_pe+onx-1
       vrho(:,j,:) = vrho(:,j,:)+v(:,j,:,tau)*0.5*( rho(:,j,:,tau)+rho(:,j+1,:,tau) )
      enddo
      do k=1,nz-1
        wrho(:,:,k) = wrho(:,:,k)+w(:,:,k,tau)*0.5*( rho(:,:,k,tau)+rho(:,:,k+1,tau) )
      enddo
end subroutine diag_variance






 subroutine diag_variance_write
!======================================================================
!=======================================================================
     use main_module   
     use diag_variance_module
     implicit none
     integer :: i,j,k,n
!  calculate means
     n=max(1,nitts)
     umean = umean*maskU/n
     vmean = vmean*maskV/n
     wmean = wmean*maskW/n
     rhom  = rhom*maskT/n
!  EKE
     uvar = uvar/n - umean**2
     vvar = vvar/n - vmean**2
     call ugrid_to_tgrid(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,uvar,uvar)  
     call vgrid_to_tgrid(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,vvar,vvar)
     uvar = 0.5*( uvar+vvar)
!  Fluxes
     do i=is_pe-onx,ie_pe+onx-1
       urho(i,:,:) = urho(i,:,:)/n - umean(i,:,:)*0.5*( rhom(i,:,:)+rhom(i+1,:,:) )
     enddo
     do j=js_pe-onx,je_pe+onx-1
       vrho(:,j,:) = vrho(:,j,:)/n - vmean(:,j,:)*0.5*( rhom(:,j,:)+rhom(:,j+1,:) )
     enddo
     do k=1,nz-1
       wrho(:,:,k) = wrho(:,:,k)/n + wmean(:,:,k)*0.5*( rhom(:,:,k)+rhom(:,:,k+1) )
     enddo
!  apply  land mask
     where (maskU==0) 
        umean=spval; urho=spval
     end where
     where (maskV==0) 
        vmean=spval; vrho=spval
     end where
     where (maskW==0) 
        wmean=spval;  wrho = spval
     end where
     where (maskT==0) 
        rhom=spval; uvar=spval
     end where

     call diag_variance_write_cdf
     nitts = 0
     umean=0;vmean=0;wmean=0;rhom=0; uvar=0;vvar=0;
     urho=0; vrho=0; wrho=0
end subroutine diag_variance_write



subroutine diag_variance_write_cdf
 use main_module   
 use diag_variance_module
 implicit none
 include "netcdf.inc"
 character (len=80) :: file
 integer :: ncid,iret,k
 integer :: tdim,ilen,timeid,varid,dims(4)
 integer :: xtdim,xudim,ytdim,yudim,ztdim,zudim
 real*8 :: time
 character :: name*64, unit*32
 real*8 :: bloc(nx,ny)

 write(file,'(a,i12,a)')  'variances_',itt,'.cdf'
 call replace_space_zero(file)
 if (my_pe==0)   print'(2a)',' writing variances to file ',file(1:len_trim(file))
 call def_grid_cdf(file)
 
 if (my_pe==0) then 
    iret=nf_open(file,NF_WRITE,ncid)
    iret=nf_set_fill(ncid, NF_NOFILL, iret)
    call ncredf(ncid, iret)
    iret=nf_inq_dimid(ncid,'xt',xtdim)
    iret=nf_inq_dimid(ncid,'xu',xudim)
    iret=nf_inq_dimid(ncid,'yt',ytdim)
    iret=nf_inq_dimid(ncid,'yu',yudim)
    iret=nf_inq_dimid(ncid,'zt',ztdim)
    iret=nf_inq_dimid(ncid,'zu',zudim)
    iret=nf_inq_dimid(ncid,'Time',tdim)

    !dims = (/xudim, ytdim,ztdim, tdim/)
    !varid = ncvdef (ncid,'um', NCFLOAT,4,dims,iret)
    !name='Mean zonal velocity'; unit = 'm/s'
    !call dvcdf(ncid,varid,name,64,unit,32,spval)

    !dims = (/xtdim, yudim,ztdim, tdim/)
    !varid = ncvdef (ncid,'vm', NCFLOAT,4,dims,iret)
    !name='Mean meridional velocity'; unit = 'm/s'
    !call dvcdf(ncid,varid,name,64,unit,32,spval)

    dims = (/xtdim, ytdim,ztdim, tdim/)
    varid = ncvdef (ncid,'rhom', NCFLOAT,4,dims,iret)
    name='Mean density'; unit = 'kg/m^3'
    call dvcdf(ncid,varid,name,64,unit,32,spval)

    dims = (/xtdim, ytdim,ztdim, tdim/)
    varid = ncvdef (ncid,'eke', NCFLOAT,4,dims,iret)
    name='Eddy kinetic energy'; unit = 'm^2/s^2'
    call dvcdf(ncid,varid,name,64,unit,32,spval)

    dims = (/xudim, ytdim,ztdim, tdim/)
    varid = ncvdef (ncid,'urho', NCFLOAT,4,dims,iret)
    name='Zonal eddy density flux'; unit = 'kg/m^2/s'
    call dvcdf(ncid,varid,name,64,unit,32,spval)

    dims = (/xtdim, yudim,ztdim, tdim/)
    varid = ncvdef (ncid,'vrho', NCFLOAT,4,dims,iret)
    name='Meridional eddy density flux'; unit = 'kg/m^2/s'
    call dvcdf(ncid,varid,name,64,unit,32,spval)

    dims = (/xtdim, ytdim,zudim, tdim/)
    varid = ncvdef (ncid,'wrho', NCFLOAT,4,dims,iret)
    name='Vertical eddy density flux'; unit = 'kg/m^2/s'
    call dvcdf(ncid,varid,name,64,unit,32,spval)

    iret= nf_put_att_int(ncid,nf_global,'nitts',nf_int,1,nitts)
    call ncendf(ncid, iret)
    call ncclos (ncid, iret)
 
    iret=nf_open(file,NF_WRITE,ncid)
    iret=nf_set_fill(ncid, NF_NOFILL, iret)
    iret=nf_inq_dimid(ncid,'Time',tdim)
    iret=nf_inq_dimlen(ncid,tdim,ilen)
    iret=nf_inq_varid(ncid,'Time',timeid)
    ilen=ilen+1
    time = itt*dt_tracer/86400 
    iret= nf_put_vara_double(ncid,timeid,(/ilen/),(/1/),(/time/))
 endif ! my_pe ==0

 do k=1,nz 
      !bloc(is_pe:ie_pe,js_pe:je_pe) = umean(is_pe:ie_pe,js_pe:je_pe,k)
      !call pe0_recv_2D(nx,ny,bloc)
      !if (my_pe==0) then
      !  iret=nf_inq_varid(ncid,'um',varid)
      !  iret= nf_put_vara_double(ncid,varid,(/1,1,k,1/), (/nx,ny,1,1/),bloc)
      !endif ! my_pe ==0

      !bloc(is_pe:ie_pe,js_pe:je_pe) = vmean(is_pe:ie_pe,js_pe:je_pe,k)
      !call pe0_recv_2D(nx,ny,bloc)
      !if (my_pe==0) then
      !  iret=nf_inq_varid(ncid,'vm',varid)
      !  iret= nf_put_vara_double(ncid,varid,(/1,1,k,1/), (/nx,ny,1,1/),bloc)
      !endif ! my_pe ==0

      bloc(is_pe:ie_pe,js_pe:je_pe) = rhom(is_pe:ie_pe,js_pe:je_pe,k)
      call pe0_recv_2D(nx,ny,bloc)
      if (my_pe==0) then
        iret=nf_inq_varid(ncid,'rhom',varid)
        iret= nf_put_vara_double(ncid,varid,(/1,1,k,1/), (/nx,ny,1,1/),bloc)
      endif ! my_pe ==0

      bloc(is_pe:ie_pe,js_pe:je_pe) = uvar(is_pe:ie_pe,js_pe:je_pe,k)
      call pe0_recv_2D(nx,ny,bloc)
      if (my_pe==0) then
        iret=nf_inq_varid(ncid,'eke',varid)
        iret= nf_put_vara_double(ncid,varid,(/1,1,k,1/), (/nx,ny,1,1/),bloc)
      endif ! my_pe ==0

      bloc(is_pe:ie_pe,js_pe:je_pe) = urho(is_pe:ie_pe,js_pe:je_pe,k)
      call pe0_recv_2D(nx,ny,bloc)
      if (my_pe==0) then
        iret=nf_inq_varid(ncid,'urho',varid)
        iret= nf_put_vara_double(ncid,varid,(/1,1,k,1/), (/nx,ny,1,1/),bloc)
      endif ! my_pe ==0

      bloc(is_pe:ie_pe,js_pe:je_pe) = vrho(is_pe:ie_pe,js_pe:je_pe,k)
      call pe0_recv_2D(nx,ny,bloc)
      if (my_pe==0) then
        iret=nf_inq_varid(ncid,'vrho',varid)
        iret= nf_put_vara_double(ncid,varid,(/1,1,k,1/), (/nx,ny,1,1/),bloc)
      endif ! my_pe ==0

      bloc(is_pe:ie_pe,js_pe:je_pe) = wrho(is_pe:ie_pe,js_pe:je_pe,k)
      call pe0_recv_2D(nx,ny,bloc)
      if (my_pe==0) then
        iret=nf_inq_varid(ncid,'wrho',varid)
        iret= nf_put_vara_double(ncid,varid,(/1,1,k,1/), (/nx,ny,1,1/),bloc)
      endif ! my_pe ==0

 enddo
         
 if (my_pe==0)  call ncclos (ncid, iret)
end subroutine diag_variance_write_cdf









subroutine diag_variance_read_restart
!=======================================================================
! read unfinished variances from file
!=======================================================================
 use main_module
 use diag_variance_module
 implicit none
 character (len=80) :: filename
 logical :: file_exists
 integer :: io,nx_,ny_,nz_,is_,ie_,js_,je_,is,ie,js,je,ierr

 is=is_pe-onx; ie=ie_pe+onx; js=js_pe-onx; je=je_pe+onx
 write(filename,'(a,i5,a)')  'unfinished_variance_PE_',my_pe,'.dta'
 call replace_space_zero(filename)
 inquire ( FILE=filename, EXIST=file_exists )
 if (.not. file_exists) then
       if (my_pe==0) then
         print'(a,a,a)',' file ',filename(1:len_trim(filename)),' not present'
         print'(a)',' reading no unfinished variances'
       endif
       return
 endif

 if (my_pe==0) print'(2a)',' reading unfinished variances from ',filename(1:len_trim(filename))
 call get_free_iounit(io,ierr)
 if (ierr/=0) goto 10
 open(io,file=filename,form='unformatted',status='old',err=10)
 read(io,err=10) nx_,ny_,nz_
 if (nx/=nx_ .or. ny/=ny_ .or. nz/= nz_) then 
       if (my_pe==0) then
        print*,' read dimensions: ',nx_,ny_,nz_
        print*,' does not match dimensions   : ',nx,ny,nz
       endif
       goto 10
 endif
 read(io,err=10) is_,ie_,js_,je_
 if (is_/=is.or.ie_/=ie.or.js_/=js.or.je_/=je) then
       if (my_pe==0) then
        print*,' read PE boundaries   ',is_,ie_,js_,je_
        print*,' which does not match ',is,ie,js,je
       endif
       goto 10
 endif
 read(io,err=10) nitts
 read(io,err=10) umean(:,:,:)
 read(io,err=10) vmean(:,:,:)
 read(io,err=10) wmean(:,:,:)
 read(io,err=10) rhom(:,:,:)
 read(io,err=10) uvar(:,:,:)
 read(io,err=10) vvar(:,:,:)
 read(io,err=10) urho(:,:,:)
 read(io,err=10) vrho(:,:,:)
 read(io,err=10) wrho(:,:,:)
 close(io)
 call fortran_barrier()
 return
 10 continue
 print'(a)',' Warning: error reading file'
end subroutine diag_variance_read_restart





subroutine diag_variance_write_restart
!=======================================================================
! write unfinished variances to restart file
!=======================================================================
 use main_module
 use diag_variance_module
 implicit none
 character (len=80) :: filename
 integer :: io,is,ie,js,je,ierr

 is=is_pe-onx; ie=ie_pe+onx; js=js_pe-onx; je=je_pe+onx
 write(filename,'(a,i5,a)')  'unfinished_variance_PE_',my_pe,'.dta'
 call replace_space_zero(filename)

 if (my_pe==0) print'(2a)',' writing unfinished variances to ',filename(1:len_trim(filename))
 call get_free_iounit(io,ierr)
 if (ierr/=0) goto 10
 open(io,file=filename,form='unformatted',status='unknown')
 write(io,err=10) nx,ny,nz
 write(io,err=10) is,ie,js,je
 write(io,err=10) nitts
 write(io,err=10) umean(is:ie,js:je,:)
 write(io,err=10) vmean(is:ie,js:je,:)
 write(io,err=10) wmean(is:ie,js:je,:)
 write(io,err=10) rhom(is:ie,js:je,:)
 write(io,err=10) uvar(is:ie,js:je,:)
 write(io,err=10) vvar(is:ie,js:je,:)
 write(io,err=10) urho(is:ie,js:je,:)
 write(io,err=10) vrho(is:ie,js:je,:)
 write(io,err=10) wrho(is:ie,js:je,:)
 close(io)
 call fortran_barrier()
 return
 10 continue
 print'(a)',' Warning: error writing file'
end subroutine diag_variance_write_restart


