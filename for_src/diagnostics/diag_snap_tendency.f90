
subroutine diag_snap_tendency
 use main_module
 implicit none
 include "netcdf.inc"
 integer :: ncid,iret,itdimid,ilen,itimeid,id,k
 real*8 :: cloc(nx,ny),fxa
 real*8, parameter :: spval = -1.0d33

 if (my_pe == 0 ) then
   iret=nf_open('pyOM_tendency.cdf',NF_WRITE,ncid)
   print'(a)',' writing to file pyOM_tendency.cdf'
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
    iret=nf_open('pyOM_tendency.cdf',NF_WRITE,ncid)
    iret=nf_inq_dimid(ncid,'Time',itdimid)
    iret=nf_inq_dimlen(ncid, itdimid,ilen)
 endif

 do k=1,nz  
    cloc(is_pe:ie_pe,js_pe:je_pe) = dtemp(is_pe:ie_pe,js_pe:je_pe,k,tau)
    where( maskT(is_pe:ie_pe,js_pe:je_pe,k) == 0.) cloc(is_pe:ie_pe,js_pe:je_pe) = spval
    call pe0_recv_2D(nx,ny,cloc)
    if (my_pe == 0 ) then      
         iret=nf_inq_varid(ncid,'dtemp_adv',id)
         iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),cloc)
    endif
 
    cloc(is_pe:ie_pe,js_pe:je_pe) = dsalt(is_pe:ie_pe,js_pe:je_pe,k,tau)
    where( maskT(is_pe:ie_pe,js_pe:je_pe,k) == 0.) cloc(is_pe:ie_pe,js_pe:je_pe) = spval
    call pe0_recv_2D(nx,ny,cloc)
    if (my_pe == 0 ) then      
         iret=nf_inq_varid(ncid,'dsalt_adv',id)
         iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),cloc)
    endif

    cloc(is_pe:ie_pe,js_pe:je_pe) = dtemp_vmix(is_pe:ie_pe,js_pe:je_pe,k)
    where( maskT(is_pe:ie_pe,js_pe:je_pe,k) == 0.) cloc(is_pe:ie_pe,js_pe:je_pe) = spval
    call pe0_recv_2D(nx,ny,cloc)
    if (my_pe == 0 ) then      
         iret=nf_inq_varid(ncid,'dtemp_vmix',id)
         iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),cloc)
    endif
 
    cloc(is_pe:ie_pe,js_pe:je_pe) = dsalt_vmix(is_pe:ie_pe,js_pe:je_pe,k)
    where( maskT(is_pe:ie_pe,js_pe:je_pe,k) == 0.) cloc(is_pe:ie_pe,js_pe:je_pe) = spval
    call pe0_recv_2D(nx,ny,cloc)
    if (my_pe == 0 ) then      
         iret=nf_inq_varid(ncid,'dsalt_vmix',id)
         iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),cloc)
    endif
 
    cloc(is_pe:ie_pe,js_pe:je_pe) = dtemp_iso(is_pe:ie_pe,js_pe:je_pe,k)
    where( maskT(is_pe:ie_pe,js_pe:je_pe,k) == 0.) cloc(is_pe:ie_pe,js_pe:je_pe) = spval
    call pe0_recv_2D(nx,ny,cloc)
    if (my_pe == 0 ) then      
         iret=nf_inq_varid(ncid,'dtemp_iso',id)
         iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),cloc)
    endif
 
    cloc(is_pe:ie_pe,js_pe:je_pe) = dsalt_iso(is_pe:ie_pe,js_pe:je_pe,k)
    where( maskT(is_pe:ie_pe,js_pe:je_pe,k) == 0.) cloc(is_pe:ie_pe,js_pe:je_pe) = spval
    call pe0_recv_2D(nx,ny,cloc)
    if (my_pe == 0 ) then      
         iret=nf_inq_varid(ncid,'dsalt_iso',id)
         iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),cloc)
    endif
 
 enddo

 if (my_pe==0)   iret = nf_close (ncid)

end subroutine diag_snap_tendency


subroutine init_diag_snap_tendency
 use main_module
 include "netcdf.inc"
 integer :: ncid,iret,lon_tdim,lon_udim,itimedim
 integer :: lat_tdim,lat_udim,id,z_tdim,z_udim, dims(4)
 character :: name*60, unit*16
 real*8, parameter :: spval = -1.0d33

 if (my_pe==0) print'(a)',' preparing file pyOM_tendency.cdf'

 call def_grid_cdf('pyOM_tendency.cdf')

 if (my_pe==0) then
      iret=nf_open('pyOM_tendency.cdf',NF_WRITE, ncid)
      iret=nf_set_fill(ncid, NF_NOFILL, iret)
      call ncredf(ncid, iret)
      iret=nf_inq_dimid(ncid,'xt',lon_tdim)
      iret=nf_inq_dimid(ncid,'xu',lon_udim)
      iret=nf_inq_dimid(ncid,'yt',lat_tdim)
      iret=nf_inq_dimid(ncid,'yu',lat_udim)
      iret=nf_inq_dimid(ncid,'zt',z_tdim)
      iret=nf_inq_dimid(ncid,'zu',z_udim)
      iret=nf_inq_dimid(ncid,'Time',itimedim)

      dims = (/Lon_tdim,lat_tdim,z_tdim,iTimedim/)
      id  = ncvdef (ncid,'dtemp_adv', NCFLOAT,4,dims,iret)
      name = 'Temperature tendency by advection'; unit = 'K/s'
      call dvcdf(ncid,id,name,60,unit,16,spval)

      dims = (/Lon_tdim,lat_tdim,z_tdim,iTimedim/)
      id  = ncvdef (ncid,'dsalt_adv', NCFLOAT,4,dims,iret)
      name = 'Salinity tendency by advection'; unit = 'g/Kg/s'
      call dvcdf(ncid,id,name,60,unit,16,spval)
      
      dims = (/Lon_tdim,lat_tdim,z_tdim,iTimedim/)
      id  = ncvdef (ncid,'dtemp_vmix', NCFLOAT,4,dims,iret)
      name = 'Temperature tendency by vert. mixing'; unit = 'K/s'
      call dvcdf(ncid,id,name,60,unit,16,spval)

      dims = (/Lon_tdim,lat_tdim,z_tdim,iTimedim/)
      id  = ncvdef (ncid,'dsalt_vmix', NCFLOAT,4,dims,iret)
      name = 'Salinity tendency by vert. mixing'; unit = 'g/Kg/s'
      call dvcdf(ncid,id,name,60,unit,16,spval)
      
      dims = (/Lon_tdim,lat_tdim,z_tdim,iTimedim/)
      id  = ncvdef (ncid,'dtemp_iso', NCFLOAT,4,dims,iret)
      name = 'Temperature tendency by isop. mixing'; unit = 'K/s'
      call dvcdf(ncid,id,name,60,unit,16,spval)

      dims = (/Lon_tdim,lat_tdim,z_tdim,iTimedim/)
      id  = ncvdef (ncid,'dsalt_iso', NCFLOAT,4,dims,iret)
      name = 'Salinity tendency by isop. mixing'; unit = 'g/Kg/s'
      call dvcdf(ncid,id,name,60,unit,16,spval)
      
         
      call ncclos (ncid, iret)
      print'(a)',' done '
 endif
 call fortran_barrier


end subroutine init_diag_snap_tendency