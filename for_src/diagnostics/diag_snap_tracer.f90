


subroutine init_snap_cdf_tracer
!=======================================================================
!     initialize NetCDF snapshot file for tracer variables
!=======================================================================
 use main_module   
 use tracer_module   
 use diagnostics_module
 implicit none
 include "netcdf.inc"
 integer :: ncid,iret
 integer :: lon_tdim,lon_udim,itimedim
 integer :: lat_tdim,lat_udim,id,z_tdim,z_udim
 integer :: dims(4),n
 character :: name*32, unit*16
 real*8, parameter :: spval = -1.0d33

 if (my_pe==0) print'(a)',' preparing file pyOM_tracer.cdf'

 call def_grid_cdf('pyOM_tracer.cdf')

 if (my_pe==0) then
      iret=nf_open('pyOM_tracer.cdf',NF_WRITE, ncid)
      iret=nf_set_fill(ncid, NF_NOFILL, iret)
      call ncredf(ncid, iret)
      iret=nf_inq_dimid(ncid,'xt',lon_tdim)
      iret=nf_inq_dimid(ncid,'xu',lon_udim)
      iret=nf_inq_dimid(ncid,'yt',lat_tdim)
      iret=nf_inq_dimid(ncid,'yu',lat_udim)
      iret=nf_inq_dimid(ncid,'zt',z_tdim)
      iret=nf_inq_dimid(ncid,'zu',z_udim)
      iret=nf_inq_dimid(ncid,'Time',itimedim)

      do n=1,n_trac
       write(name, '("trac",i3)') n; call replace_space_zero(name)
       dims = (/Lon_tdim,lat_tdim,z_udim,iTimedim/)
       id  = ncvdef (ncid,name, NCFLOAT,4,dims,iret)
       name = 'passive tracer'; unit = ' '
       call dvcdf(ncid,id,name,32,unit,16,spval)

       write(name, '("trac",i3,"_sflux")') n; call replace_space_zero(name)
       dims = (/Lon_tdim,lat_tdim,iTimedim,1/)
       id  = ncvdef (ncid,name, NCFLOAT,3,dims,iret)
       name = 'surface flux'; unit = ' '
       call dvcdf(ncid,id,name,len_trim(name),unit,16,spval)
      enddo
      call ncclos (ncid, iret)
      print'(a)',' done '
 endif
 call fortran_barrier
end subroutine init_snap_cdf_tracer


subroutine diag_snap_tracer
!=======================================================================
!     write tracer variables to NetCDF snapshot file
!=======================================================================
 use main_module   
 use tracer_module   
 use diagnostics_module
 implicit none
 include "netcdf.inc"
 integer :: ncid,iret,n
 real*8 :: bloc(nx,ny),fxa
 integer :: itdimid,ilen,itimeid,id,k
 real*8, parameter :: spval = -1.0d33
 character :: name*32

 if (my_pe == 0 ) then
   iret=nf_open('pyOM_tracer.cdf',NF_WRITE,ncid)
   print'(a)',' writing to file pyOM_tracer.cdf'
   iret=nf_set_fill(ncid, NF_NOFILL, iret)
   iret=nf_inq_dimid(ncid,'Time',itdimid)
   iret=nf_inq_dimlen(ncid, itdimid,ilen)
   ilen=ilen+1
   fxa = itt*dt_tracer/86400.0
   iret=nf_inq_varid(ncid,'Time',itimeid)
   iret= nf_put_vara_double(ncid,itimeid,ilen,1,fxa)
   call ncclos (ncid, iret)
 endif
 call fortran_barrier

 if (my_pe == 0 ) then
    iret=nf_open('pyOM_tracer.cdf',NF_WRITE,ncid)
    iret=nf_inq_dimid(ncid,'Time',itdimid)
    iret=nf_inq_dimlen(ncid, itdimid,ilen)
 endif

 do n=1,n_trac
   bloc(is_pe:ie_pe,js_pe:je_pe) = forc_trac_surface(is_pe:ie_pe,js_pe:je_pe,n)
   where( maskW(is_pe:ie_pe,js_pe:je_pe,nz) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
   call pe0_recv_2D(nx,ny,bloc)
   if (my_pe==0) then
     write(name, '("trac",i3,"_sflux")') n; call replace_space_zero(name)
    iret=nf_inq_varid(ncid,name,id)
    iret= nf_put_vara_double(ncid,id,(/1,1,ilen/), (/nx,ny,1/),bloc)
   endif
 enddo
 ! 3D fields

 do k=1,nz
   do n=1,n_trac
      bloc(is_pe:ie_pe,js_pe:je_pe) = trac(is_pe:ie_pe,js_pe:je_pe,k,tau,n)
      where( maskW(is_pe:ie_pe,js_pe:je_pe,k) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
      call pe0_recv_2D(nx,ny,bloc)
      if (my_pe == 0 ) then 
         write(name, '("trac",i3)') n; call replace_space_zero(name)
         iret=nf_inq_varid(ncid,name,id)
         iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),bloc)
      endif
   enddo
 enddo

 if (my_pe==0)   iret = nf_close (ncid)
end subroutine diag_snap_tracer





